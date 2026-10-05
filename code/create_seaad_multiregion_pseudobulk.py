#!/usr/bin/env python3
"""Create a library-by-Supertype SEA-AD Multiregion pseudobulk H5AD.

The script uses the public ABC Atlas Access manifest and the ten regional raw
count matrices in ``SEA-AD-Multiregion-10X``. It processes one regional H5AD at
a time in backed mode and combines the regional pseudobulk results.

Output
------
One H5AD with:

* one observation per ``library_label x Supertype``;
* summed integer raw counts in ``X``;
* expressing-nucleus counts in ``layers['number_of_nuclei_expressing']``;
* total contributing nuclei in ``obs['number_of_nuclei']``;
* cell taxonomy, library, donor, and disease metadata in ``obs``;
* source gene metadata in ``var``.

Nuclei without a valid Supertype assignment are reported and excluded. The
source raw H5AD files are downloaded/resolved by ``AbcProjectCache`` and opened
in backed mode; they are never loaded in full into RAM.

Suggested Code Ocean environment
--------------------------------
pip install "abc_atlas_access[notebooks] @ git+https://github.com/AllenInstitute/abc_atlas_access.git" \
            anndata pandas numpy scipy h5py

Example
-------
python code/create_seaad_multiregion_pseudobulk.py \
  --cache-dir /tmp/abc_atlas_cache \
  --output /results/SEA_AD_Multiregion_library_supertype_pseudobulk.h5ad
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache


DATASET_DIR = "SEA-AD-Multiregion-10X"
TAXONOMY_DIR = "SEA-AD-Multiregion-taxonomy"
DEFAULT_MANIFEST = "releases/20260711/manifest.json"
DEFAULT_OUTPUT = "/results/SEA_AD_Multiregion_library_supertype_pseudobulk.h5ad"
DEFAULT_REGIONS = (
    "AnG",
    "DFC",
    "FI",
    "HIP",
    "ITG",
    "LEC",
    "MEC",
    "MTG",
    "STG",
    "V1C",
)
GROUPING_COLUMNS = ("library_label", "Supertype")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create summed-count and expressing-nucleus pseudobulks for the "
            "SEA-AD Multiregion dataset, grouped by library and Supertype."
        )
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("/tmp/abc_atlas_cache"),
        help=(
            "ABC Atlas download cache. Regional raw H5ADs are resolved through "
            "this cache (default: /tmp/abc_atlas_cache)."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT),
        help=f"Output H5AD (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--manifest",
        default=DEFAULT_MANIFEST,
        help=(
            "ABC Atlas manifest. Use 'latest' to retain the cache default "
            f"(default: {DEFAULT_MANIFEST})."
        ),
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        default=list(DEFAULT_REGIONS),
        help="Region prefixes to process (default: all ten SEA-AD regions).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=4096,
        help="Source nuclei read per matrix chunk (default: 4096).",
    )
    return parser.parse_args()


def load_manifest(cache: AbcProjectCache, requested: str) -> None:
    """Load the requested ABC Atlas manifest."""
    if requested.lower() == "latest":
        logging.info("Using cache default manifest: %s", cache.current_manifest)
        return

    if requested not in cache.cache.manifest_file_names:
        cache.cache.manifest_file_names.append(requested)
    cache.load_manifest(requested)
    logging.info("Using manifest: %s", cache.current_manifest)


def load_metadata(cache: AbcProjectCache) -> pd.DataFrame:
    """Build cell metadata with taxonomy, library, donor, and disease fields."""
    logging.info("Loading cell metadata")
    cell = cache.get_metadata_dataframe(
        directory=DATASET_DIR,
        file_name="cell_metadata",
        low_memory=False,
    ).set_index("cell_label", drop=True)
    cell.index = cell.index.astype(str)
    cell.index.name = "cell_label"
    if not cell.index.is_unique:
        raise ValueError("cell_metadata contains duplicate cell_label values")

    logging.info("Loading library, donor, and disease metadata")
    library = cache.get_metadata_dataframe(
        directory=DATASET_DIR,
        file_name="library",
        low_memory=False,
    ).set_index("library_label", drop=True)
    donor = cache.get_metadata_dataframe(
        directory=DATASET_DIR,
        file_name="donor",
        low_memory=False,
    ).set_index("donor_label", drop=True)
    disease = cache.get_metadata_dataframe(
        directory=DATASET_DIR,
        file_name="disease",
        low_memory=False,
    ).set_index("donor_label", drop=True)

    library.index = library.index.astype(str)
    donor.index = donor.index.astype(str)
    disease.index = disease.index.astype(str)

    logging.info("Loading taxonomy metadata")
    cell_membership = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIR,
        file_name="cell_to_cluster_membership",
        low_memory=False,
    ).set_index("cell_label", drop=True)
    cell_membership.index = cell_membership.index.astype(str)
    if "cluster_alias" not in cell_membership.columns:
        raise KeyError("cell_to_cluster_membership lacks cluster_alias")
    if not cell_membership.index.is_unique:
        raise ValueError(
            "cell_to_cluster_membership contains duplicate cell_label values"
        )

    cluster = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIR,
        file_name="cluster",
        dtype={"number_of_cells": "Int64"},
    )
    required_cluster = {"label", "cluster_alias"}
    missing_cluster = required_cluster.difference(cluster.columns)
    if missing_cluster:
        raise KeyError(
            "cluster.csv lacks required columns: "
            + ", ".join(sorted(missing_cluster))
        )

    term_sets = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIR,
        file_name="cluster_annotation_term_set",
    )
    if not {"name", "order"}.issubset(term_sets.columns):
        raise KeyError("cluster_annotation_term_set lacks name or order")
    term_sets = term_sets.sort_values("order", kind="stable")
    taxonomy_levels = term_sets["name"].astype(str).tolist()
    expected_levels = ["Class", "Subclass", "Supertype"]
    if taxonomy_levels != expected_levels:
        raise ValueError(
            "Unexpected SEA-AD taxonomy levels. Expected "
            f"{expected_levels!r}; found {taxonomy_levels!r}"
        )

    terms = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIR,
        file_name="cluster_annotation_term",
    ).rename(columns={"label": "cluster_annotation_term_label"})
    required_terms = {
        "cluster_annotation_term_label",
        "cluster_annotation_term_set_name",
        "name",
        "color_hex_triplet",
        "term_order",
    }
    missing_terms = required_terms.difference(terms.columns)
    if missing_terms:
        raise KeyError(
            "cluster_annotation_term lacks required columns: "
            + ", ".join(sorted(missing_terms))
        )
    terms = terms.set_index("cluster_annotation_term_label", drop=True)

    cluster_membership = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIR,
        file_name="cluster_to_cluster_annotation_membership",
        low_memory=False,
    ).set_index("cluster_annotation_term_label", drop=True)
    if "cluster_alias" not in cluster_membership.columns:
        raise KeyError(
            "cluster_to_cluster_annotation_membership lacks cluster_alias"
        )

    membership_with_terms = cluster_membership.join(
        terms,
        how="left",
        rsuffix="_annotation_term",
    ).reset_index()

    required_expanded = {
        "cluster_alias",
        "cluster_annotation_term_set_name",
        "name",
        "color_hex_triplet",
        "term_order",
    }
    missing_expanded = required_expanded.difference(membership_with_terms.columns)
    if missing_expanded:
        raise KeyError(
            "Expanded taxonomy membership lacks: "
            + ", ".join(sorted(missing_expanded))
        )

    grouped = membership_with_terms.groupby(
        ["cluster_alias", "cluster_annotation_term_set_name"],
        observed=True,
        sort=False,
    )
    cluster_details = grouped["name"].first().unstack()
    cluster_colors = grouped["color_hex_triplet"].first().unstack()
    cluster_orders = grouped["term_order"].first().unstack()

    for frame, suffix in (
        (cluster_details, ""),
        (cluster_colors, "_color"),
        (cluster_orders, "_order"),
    ):
        missing_levels = set(expected_levels).difference(frame.columns)
        if missing_levels:
            raise ValueError(
                f"Taxonomy lookup {suffix!r} lacks levels: "
                + ", ".join(sorted(missing_levels))
            )

    cluster_details = cluster_details[expected_levels]
    cluster_colors = cluster_colors[expected_levels].add_suffix("_color")
    cluster_orders = cluster_orders[expected_levels].add_suffix("_order")
    cluster_lookup = cluster_details.join(cluster_colors).join(cluster_orders)

    # Reproduce the joins in the published SEA-AD notebooks.
    extended = cell.join(
        cell_membership[["cluster_alias"]],
        how="inner",
    )
    extended = extended.join(cluster_lookup, on="cluster_alias", how="left")
    extended = extended.join(
        library,
        on="library_label",
        how="left",
        rsuffix="_library",
    )

    donor_join_column = (
        "donor_label_library"
        if "donor_label_library" in extended.columns
        else "donor_label"
    )
    if donor_join_column not in extended.columns:
        raise KeyError("Expanded metadata lacks a donor label for donor joins")

    extended = extended.join(
        donor,
        on=donor_join_column,
        how="left",
        rsuffix="_donor",
    )
    extended = extended.join(
        disease,
        on=donor_join_column,
        how="left",
        rsuffix="_disease",
    )

    required = {
        "library_label",
        "cluster_alias",
        "Class",
        "Subclass",
        "Supertype",
    }
    missing = required.difference(extended.columns)
    if missing:
        raise KeyError(
            "Expanded cell metadata lacks required columns: "
            + ", ".join(sorted(missing))
        )

    for column in required:
        if extended[column].isna().any():
            raise ValueError(f"Expanded metadata contains missing {column} values")

    extended["library_label"] = extended["library_label"].astype(str)
    extended["cluster_alias"] = extended["cluster_alias"].astype(str)
    extended["Supertype"] = extended["Supertype"].astype(str)

    validate_library_invariance(extended)

    logging.info(
        "Expanded metadata: %s taxonomy-assigned nuclei, %s libraries, "
        "%s supertypes",
        len(extended),
        extended["library_label"].nunique(),
        extended["Supertype"].nunique(),
    )
    return extended


def validate_library_invariance(metadata: pd.DataFrame) -> None:
    """Ensure each library maps to one donor and one brain region."""
    candidate_columns = [
        column
        for column in (
            "donor_label",
            "donor_label_library",
            "Brain Region",
            "region_of_interest_label",
            "feature_matrix_label",
        )
        if column in metadata.columns
    ]
    for column in candidate_columns:
        number_per_library = metadata.groupby(
            "library_label",
            observed=True,
        )[column].nunique(dropna=False)
        bad = number_per_library[number_per_library > 1]
        if not bad.empty:
            raise ValueError(
                f"{len(bad)} libraries map to multiple {column} values; "
                f"examples: {bad.index[:10].tolist()}"
            )


def expression_matrix_names(
    cache: AbcProjectCache,
    regions: Sequence[str],
) -> list[str]:
    """Resolve one raw-count expression entry for every requested region."""
    if hasattr(cache, "list_expression_matrix_files"):
        available: Iterable[str] = cache.list_expression_matrix_files(
            directory=DATASET_DIR
        )
    else:
        available = cache.list_data_files(directory=DATASET_DIR)
    available = [str(name).strip("/") for name in available]

    result: list[str] = []
    for region in regions:
        expected = f"{region}-10X/raw"
        matches = [name for name in available if name == expected]
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected exactly one expression entry {expected!r}; "
                f"found {matches!r}. Available entries: {available!r}"
            )
        result.append(matches[0])
    return result


def gene_metadata_for_source(
    cache: AbcProjectCache,
    source: ad.AnnData,
) -> pd.DataFrame:
    """Create var metadata in source order, indexed by stable unique gene key."""
    gene = cache.get_metadata_dataframe(
        directory=DATASET_DIR,
        file_name="gene",
        low_memory=False,
    )
    required = {"gene_identifier", "gene_symbol"}
    missing = required.difference(gene.columns)
    if missing:
        raise KeyError("gene.csv lacks: " + ", ".join(sorted(missing)))

    gene = gene.drop_duplicates("gene_identifier").copy()
    gene["gene_identifier"] = gene["gene_identifier"].astype(str)
    gene = gene.set_index("gene_identifier", drop=False)

    source_ids = pd.Index(
        source.var_names.astype(str),
        name="gene_identifier",
    )
    missing_ids = source_ids.difference(gene.index)
    if len(missing_ids):
        raise ValueError(
            f"{len(missing_ids)} source genes are absent from gene.csv; "
            f"examples: {missing_ids[:10].tolist()}"
        )

    var = gene.reindex(source_ids).copy()
    symbols = var["gene_symbol"].astype("string")
    fallback = pd.Series(source_ids.to_numpy(), index=source_ids, dtype="string")
    symbols = symbols.fillna(fallback).astype(str)

    symbol_counts = symbols.value_counts()
    unique_keys = [
        symbol if symbol_counts[symbol] == 1 else f"{symbol}__{gene_id}"
        for symbol, gene_id in zip(symbols, source_ids)
    ]
    var.insert(0, "gene_symbol_original", symbols.to_numpy())
    var.index = pd.Index(unique_keys, name="gene_key")
    if not var.index.is_unique:
        raise ValueError("Constructed gene_key values are not unique")
    return var


def make_group_metadata(
    valid_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Create deterministic library-by-Supertype groups and metadata."""
    group_frame = valid_metadata[list(GROUPING_COLUMNS)].copy()
    group_index = pd.MultiIndex.from_frame(group_frame)
    codes, unique_groups = pd.factorize(group_index, sort=True)

    group_keys = pd.DataFrame(
        unique_groups.tolist(),
        columns=list(GROUPING_COLUMNS),
    )
    group_keys["_group_code"] = np.arange(
        len(group_keys),
        dtype=np.int64,
    )

    # Metadata columns besides taxonomy colors/orders must be invariant within a
    # library x Supertype row. Keep the first value after validating the fields
    # that define the biological and technical unit.
    required_invariant = [
        column
        for column in (
            "cluster_alias",
            "Class",
            "Subclass",
            "Supertype",
            "donor_label",
            "donor_label_library",
            "Brain Region",
            "region_of_interest_label",
            "feature_matrix_label",
        )
        if column in valid_metadata.columns
    ]
    temporary = valid_metadata.assign(_group_code=codes)
    for column in required_invariant:
        counts = temporary.groupby(
            "_group_code",
            observed=True,
        )[column].nunique(dropna=False)
        bad = counts[counts > 1]
        if not bad.empty:
            raise ValueError(
                f"Pseudobulk groups contain multiple {column} values; "
                f"group codes: {bad.index[:10].tolist()}"
            )

    first_metadata = temporary.groupby(
        "_group_code",
        sort=True,
        observed=True,
    ).first()
    first_metadata = first_metadata.reindex(group_keys["_group_code"])
    first_metadata = first_metadata.drop(
        columns=list(GROUPING_COLUMNS),
        errors="ignore",
    )

    obs = pd.concat(
        [group_keys.set_index("_group_code"), first_metadata],
        axis=1,
    )
    number_of_nuclei = np.bincount(
        codes,
        minlength=len(obs),
    ).astype(np.int64)
    obs.insert(2, "number_of_nuclei", number_of_nuclei)

    obs.index = pd.Index(
        [
            f"{library_label}__supertype_{cluster_alias}"
            for library_label, cluster_alias in zip(
                obs["library_label"],
                obs["cluster_alias"],
            )
        ],
        name="pseudobulk_id",
    )
    if not obs.index.is_unique:
        raise ValueError("Pseudobulk identifiers are not unique")
    return obs, codes.astype(np.int64, copy=False)


def aggregate_counts(
    source: ad.AnnData,
    full_group_codes: np.ndarray,
    n_groups: int,
    chunk_size: int,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Aggregate raw counts and expressing-nucleus counts in row chunks."""
    n_nuclei, n_genes = source.shape
    if len(full_group_codes) != n_nuclei:
        raise ValueError("Group-code vector does not match source row count")

    total_counts = sparse.csr_matrix((n_groups, n_genes), dtype=np.int64)
    total_detected = sparse.csr_matrix((n_groups, n_genes), dtype=np.int64)

    for start in range(0, n_nuclei, chunk_size):
        stop = min(start + chunk_size, n_nuclei)
        local_codes = full_group_codes[start:stop]
        keep = local_codes >= 0
        if not np.any(keep):
            continue

        logging.info(
            "Aggregating source rows %s-%s of %s (%s assigned nuclei)",
            start + 1,
            stop,
            n_nuclei,
            int(keep.sum()),
        )
        block = source.X[start:stop]
        if sparse.issparse(block):
            block = block.tocsr()
        else:
            block = sparse.csr_matrix(np.asarray(block))
        block = block[keep, :]

        if block.nnz and not np.allclose(block.data, np.rint(block.data)):
            raise ValueError(
                "The source matrix is not integer-valued. Confirm that the "
                "regional */raw matrix was selected instead of */log2."
            )
        block.data = np.rint(block.data).astype(np.int64, copy=False)
        if block.nnz and np.any(block.data < 0):
            raise ValueError("Raw expression matrix contains negative counts")

        detected = block.copy()
        detected.data = np.ones(detected.nnz, dtype=np.int64)

        selected_codes = local_codes[keep]
        assignment = sparse.csr_matrix(
            (
                np.ones(len(selected_codes), dtype=np.int64),
                (
                    selected_codes,
                    np.arange(len(selected_codes), dtype=np.int64),
                ),
            ),
            shape=(n_groups, len(selected_codes)),
        )
        total_counts = total_counts + assignment @ block
        total_detected = total_detected + assignment @ detected

    total_counts.sum_duplicates()
    total_counts.eliminate_zeros()
    total_detected.sum_duplicates()
    total_detected.eliminate_zeros()
    return total_counts, total_detected.astype(np.uint32)


def sanitize_dataframe_for_h5ad(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert mixed object columns to nullable strings for H5AD output."""
    result = frame.copy()
    for column in result.columns:
        if result[column].dtype == object:
            result[column] = result[column].astype("string")
    return result


def process_region(
    cache: AbcProjectCache,
    matrix_name: str,
    expanded_metadata: pd.DataFrame,
    expected_var: pd.DataFrame | None,
    chunk_size: int,
) -> tuple[pd.DataFrame, sparse.csr_matrix, sparse.csr_matrix, pd.DataFrame]:
    """Process one regional raw H5AD and return regional pseudobulk data."""
    logging.info("Resolving expression matrix: %s", matrix_name)
    source_path = Path(
        cache.get_data_path(
            directory=DATASET_DIR,
            file_name=matrix_name,
        )
    )
    logging.info("Opening regional raw matrix in backed mode: %s", source_path)
    source = ad.read_h5ad(source_path, backed="r")

    try:
        source_cells = pd.Index(source.obs_names.astype(str), name="cell_label")
        if not source_cells.is_unique:
            raise ValueError(f"{matrix_name} contains duplicate source cell labels")

        aligned = expanded_metadata.reindex(source_cells)
        assigned = aligned["library_label"].notna() & aligned["Supertype"].notna()
        n_excluded = int((~assigned).sum())
        logging.info(
            "%s: %s source nuclei; %s taxonomy-assigned; %s excluded",
            matrix_name,
            source.n_obs,
            int(assigned.sum()),
            n_excluded,
        )
        if not assigned.any():
            raise ValueError(f"{matrix_name} has no taxonomy-assigned nuclei")

        valid_metadata = aligned.loc[assigned].copy()
        expected_feature_matrix = matrix_name.split("/", maxsplit=1)[0]
        if "feature_matrix_label" in valid_metadata.columns:
            observed_labels = set(
                valid_metadata["feature_matrix_label"].dropna().astype(str)
            )
            if observed_labels != {expected_feature_matrix}:
                raise ValueError(
                    f"{matrix_name} expected feature_matrix_label "
                    f"{expected_feature_matrix!r}; found {sorted(observed_labels)!r}"
                )

        obs, selected_codes = make_group_metadata(valid_metadata)
        full_codes = np.full(source.n_obs, -1, dtype=np.int64)
        full_codes[np.flatnonzero(assigned.to_numpy())] = selected_codes

        var = gene_metadata_for_source(cache, source)
        if expected_var is not None:
            if not var.index.equals(expected_var.index):
                raise ValueError(
                    f"{matrix_name} gene keys or order differ from earlier regions"
                )
            left_ids = var.get("gene_identifier", pd.Series(var.index, index=var.index))
            right_ids = expected_var.get(
                "gene_identifier",
                pd.Series(expected_var.index, index=expected_var.index),
            )
            if not np.array_equal(
                left_ids.astype(str).to_numpy(),
                right_ids.astype(str).to_numpy(),
            ):
                raise ValueError(
                    f"{matrix_name} gene identifiers differ from earlier regions"
                )

        counts, detected = aggregate_counts(
            source=source,
            full_group_codes=full_codes,
            n_groups=len(obs),
            chunk_size=chunk_size,
        )

        detected_coo = detected.tocoo()
        row_totals = obs["number_of_nuclei"].to_numpy(dtype=np.int64)
        if np.any(
            detected_coo.data.astype(np.int64)
            > row_totals[detected_coo.row]
        ):
            raise RuntimeError(
                f"{matrix_name}: expressing-nucleus counts exceed total nuclei"
            )

        obs.insert(0, "source_expression_matrix", matrix_name)
        return obs, counts, detected, var
    finally:
        if getattr(source, "file", None) is not None:
            source.file.close()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be at least 1")

    regions = [str(region).strip() for region in args.regions]
    unknown = sorted(set(regions).difference(DEFAULT_REGIONS))
    if unknown:
        raise ValueError(
            "Unknown regions: " + ", ".join(unknown) + ". Valid regions: "
            + ", ".join(DEFAULT_REGIONS)
        )
    if len(regions) != len(set(regions)):
        raise ValueError("--regions contains duplicate values")

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    args.output = args.output.resolve()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    cache = AbcProjectCache.from_cache_dir(args.cache_dir)
    load_manifest(cache, args.manifest)
    expanded_metadata = load_metadata(cache)
    matrix_names = expression_matrix_names(cache, regions)

    all_obs: list[pd.DataFrame] = []
    count_blocks: list[sparse.csr_matrix] = []
    detected_blocks: list[sparse.csr_matrix] = []
    output_var: pd.DataFrame | None = None

    for matrix_name in matrix_names:
        obs, counts, detected, var = process_region(
            cache=cache,
            matrix_name=matrix_name,
            expanded_metadata=expanded_metadata,
            expected_var=output_var,
            chunk_size=args.chunk_size,
        )
        if output_var is None:
            output_var = var
        all_obs.append(obs)
        count_blocks.append(counts)
        detected_blocks.append(detected)

    if output_var is None:
        raise RuntimeError("No regions were processed")

    obs = pd.concat(all_obs, axis=0, verify_integrity=True)
    summed_counts = sparse.vstack(count_blocks, format="csr", dtype=np.int64)
    detected_counts = sparse.vstack(
        detected_blocks,
        format="csr",
        dtype=np.uint32,
    )

    if summed_counts.shape != detected_counts.shape:
        raise RuntimeError("Summed-count and detection-count shapes differ")
    if summed_counts.shape != (len(obs), len(output_var)):
        raise RuntimeError(
            "Combined matrix shape does not match obs and var dimensions"
        )

    obs = sanitize_dataframe_for_h5ad(obs)
    output_var = sanitize_dataframe_for_h5ad(output_var)

    pseudobulk = ad.AnnData(
        X=summed_counts,
        obs=obs,
        var=output_var,
        layers={
            "number_of_nuclei_expressing": detected_counts,
        },
        uns={
            "dataset": DATASET_DIR,
            "taxonomy": TAXONOMY_DIR,
            "manifest": str(cache.current_manifest),
            "aggregation": (
                "sum of raw counts by library_label and Supertype"
            ),
            "detection_aggregation": (
                "number of nuclei with raw count > 0 by library_label and "
                "Supertype"
            ),
            "detection_layer": "number_of_nuclei_expressing",
            "grouping_columns": list(GROUPING_COLUMNS),
            "count_column": "number_of_nuclei",
            "regional_expression_matrices": matrix_names,
            "excluded_nuclei_policy": (
                "source nuclei lacking library_label or Supertype are excluded"
            ),
        },
    )

    logging.info(
        "Writing %s pseudobulk rows x %s genes to %s",
        pseudobulk.n_obs,
        pseudobulk.n_vars,
        args.output,
    )
    pseudobulk.write_h5ad(args.output, compression="gzip")

    # Read-back validation catches serialization and schema problems.
    check = ad.read_h5ad(args.output, backed="r")
    try:
        if check.shape != pseudobulk.shape:
            raise RuntimeError(
                f"Read-back shape {check.shape} differs from {pseudobulk.shape}"
            )
        if "number_of_nuclei" not in check.obs.columns:
            raise RuntimeError("Read-back output lacks number_of_nuclei")
        if "number_of_nuclei_expressing" not in check.layers:
            raise RuntimeError(
                "Read-back output lacks number_of_nuclei_expressing layer"
            )
        if check.layers["number_of_nuclei_expressing"].shape != check.shape:
            raise RuntimeError("Read-back detection layer shape differs from X")
        for column in GROUPING_COLUMNS:
            if column not in check.obs.columns:
                raise RuntimeError(f"Read-back output lacks {column}")
    finally:
        check.file.close()

    logging.info(
        "Complete: %s rows x %s genes written to %s",
        pseudobulk.n_obs,
        pseudobulk.n_vars,
        args.output,
    )


if __name__ == "__main__":
    main()
