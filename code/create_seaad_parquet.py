#!/usr/bin/env python3
"""Memory-bounded SEA-AD Multiregion pseudobulk generation.

Processes one regional raw-count H5AD at a time, filters the two nucleus-level
metadata tables in CSV chunks, writes each regional pseudobulk to scratch, and
combines only the small regional pseudobulk objects at the end.
"""
from __future__ import annotations

import argparse, gc, logging, shutil
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

DATASET_DIR = "SEA-AD-Multiregion-10X"
TAXONOMY_DIR = "SEA-AD-Multiregion-taxonomy"
DEFAULT_MANIFEST = "releases/20260711/manifest.json"
REGIONS = ("AnG", "DFC", "FI", "HIP", "ITG", "LEC", "MEC", "MTG", "STG", "V1C")


def args():
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", type=Path, default=Path("/scratch/abc_atlas_cache"))
    p.add_argument("--scratch-dir", type=Path, default=Path("/scratch/seaad_pseudobulk"))
    p.add_argument("--output", type=Path, default=Path("/results/SEA_AD_Multiregion_library_supertype_pseudobulk.h5ad"))
    p.add_argument("--manifest", default=DEFAULT_MANIFEST)
    p.add_argument("--regions", nargs="+", default=list(REGIONS))
    p.add_argument("--matrix-chunk-size", type=int, default=1024)
    p.add_argument("--metadata-chunk-size", type=int, default=250000)
    p.add_argument("--keep-regional-files", action="store_true")
    p.add_argument(
        "--resume", action="store_true",
        help="Reuse valid regional pseudobulk checkpoints already in scratch.",
    )
    return p.parse_args()


def load_manifest(cache, name):
    if name.lower() == "latest":
        return
    if name not in cache.cache.manifest_file_names:
        cache.cache.manifest_file_names.append(name)
    cache.load_manifest(name)


def metadata_path(cache, directory, name):
    """Resolve/download a metadata CSV without loading it into pandas."""
    for method_name in ("get_metadata_path", "get_data_path"):
        method = getattr(cache, method_name, None)
        if method is not None:
            try:
                return Path(method(directory=directory, file_name=name))
            except Exception:
                pass
    # Last resort: trigger download, then ask the cache again.
    cache.get_metadata_dataframe(directory=directory, file_name=name, nrows=1)
    method = getattr(cache, "get_metadata_path", None)
    if method is None:
        raise RuntimeError(f"Cannot resolve metadata path for {directory}/{name}")
    return Path(method(directory=directory, file_name=name))


def read_small(cache, directory, name, index=None):
    df = cache.get_metadata_dataframe(directory=directory, file_name=name, low_memory=False)
    if index is not None:
        df[index] = df[index].astype(str)
        df = df.set_index(index, drop=True)
    return df


def read_filtered_csv(path, wanted_ids, id_column, chunk_size, columns=None):
    """Read only rows whose ID occurs in the current regional H5AD."""
    pieces = []
    wanted = wanted_ids if isinstance(wanted_ids, set) else set(wanted_ids)
    for chunk_no, chunk in enumerate(pd.read_csv(path, chunksize=chunk_size, low_memory=False, usecols=columns), 1):
        ids = chunk[id_column].astype(str)
        keep = ids.isin(wanted)
        if keep.any():
            out = chunk.loc[keep].copy()
            out[id_column] = ids.loc[keep].to_numpy()
            pieces.append(out)
        if chunk_no % 10 == 0:
            logging.info("Scanned %d metadata chunks from %s", chunk_no, path.name)
    if not pieces:
        return pd.DataFrame(columns=columns or [id_column]).set_index(id_column)
    result = pd.concat(pieces, ignore_index=True)
    result = result.drop_duplicates(id_column).set_index(id_column, drop=True)
    del pieces
    gc.collect()
    return result


def taxonomy_lookup(cache):
    """Build cluster_alias -> Class/Subclass/Supertype metadata robustly.

    The taxonomy release stores a term-set key on each annotation term. The
    human-readable set name (Class, Subclass, Supertype) is supplied by
    cluster_annotation_term_set.csv, so it must be joined explicitly rather
    than assumed to already exist on the merged membership table.
    """
    terms = read_small(cache, TAXONOMY_DIR, "cluster_annotation_term")
    membership = read_small(
        cache, TAXONOMY_DIR, "cluster_to_cluster_annotation_membership"
    )
    term_sets = read_small(cache, TAXONOMY_DIR, "cluster_annotation_term_set")

    # Normalize the annotation-term identifier used for the first merge.
    if "cluster_annotation_term_label" not in terms.columns:
        if "label" not in terms.columns:
            raise KeyError(
                "cluster_annotation_term.csv lacks label or "
                "cluster_annotation_term_label"
            )
        terms = terms.rename(columns={"label": "cluster_annotation_term_label"})
    if "cluster_annotation_term_label" not in membership.columns:
        raise KeyError(
            "cluster_to_cluster_annotation_membership.csv lacks "
            "cluster_annotation_term_label"
        )

    merged = membership.merge(
        terms,
        on="cluster_annotation_term_label",
        how="left",
        suffixes=("_membership", "_term"),
    )

    # Some releases carry the readable set name directly. Others carry only a
    # term-set label, which must be translated with the term-set table.
    direct_name_candidates = [
        "cluster_annotation_term_set_name",
        "cluster_annotation_term_set_name_term",
        "cluster_annotation_term_set_name_membership",
    ]
    direct_name = next(
        (column for column in direct_name_candidates if column in merged.columns),
        None,
    )

    if direct_name is not None:
        merged["taxonomy_level"] = merged[direct_name]
    else:
        set_key_candidates = [
            "cluster_annotation_term_set_label",
            "cluster_annotation_term_set_label_term",
            "cluster_annotation_term_set_label_membership",
        ]
        merged_set_key = next(
            (column for column in set_key_candidates if column in merged.columns),
            None,
        )
        term_set_key = next(
            (
                column
                for column in ("cluster_annotation_term_set_label", "label")
                if column in term_sets.columns
            ),
            None,
        )
        if merged_set_key is None or term_set_key is None or "name" not in term_sets.columns:
            raise KeyError(
                "Could not resolve taxonomy level names. "
                f"Merged columns: {merged.columns.tolist()}; "
                f"term-set columns: {term_sets.columns.tolist()}"
            )
        set_lookup = term_sets[[term_set_key, "name"]].drop_duplicates().rename(
            columns={term_set_key: "_term_set_key", "name": "taxonomy_level"}
        )
        merged = merged.merge(
            set_lookup,
            left_on=merged_set_key,
            right_on="_term_set_key",
            how="left",
        )

    # Resolve annotation-term attributes after the suffix-aware merge.
    def resolve(base):
        candidates = [base, f"{base}_term", f"{base}_membership"]
        column = next((value for value in candidates if value in merged.columns), None)
        if column is None:
            raise KeyError(
                f"Taxonomy merge lacks {base}; columns: {merged.columns.tolist()}"
            )
        return column

    name_column = resolve("name")
    color_column = resolve("color_hex_triplet")
    order_column = resolve("term_order")
    if "cluster_alias" not in merged.columns:
        raise KeyError("Taxonomy membership lacks cluster_alias")

    grouped = merged.groupby(
        ["cluster_alias", "taxonomy_level"],
        observed=True,
        sort=False,
    )
    names = grouped[name_column].first().unstack()
    colors = grouped[color_column].first().unstack().add_suffix("_color")
    orders = grouped[order_column].first().unstack().add_suffix("_order")

    required = ["Class", "Subclass", "Supertype"]
    missing = set(required).difference(names.columns)
    if missing:
        raise ValueError(
            f"Taxonomy lacks levels {sorted(missing)}; found {names.columns.tolist()}"
        )
    lookup = names[required].join(colors).join(orders)
    # cluster_alias is numeric in the taxonomy CSV but string-valued in the
    # cell-membership metadata after chunked loading. Normalize both join keys.
    lookup.index = lookup.index.astype(str)
    lookup.index.name = "cluster_alias"
    if not lookup.index.is_unique:
        raise ValueError("Taxonomy lookup contains duplicate cluster_alias values")
    return lookup


def static_metadata(cache):
    library = read_small(cache, DATASET_DIR, "library", "library_label")
    donor = read_small(cache, DATASET_DIR, "donor", "donor_label")
    disease = read_small(cache, DATASET_DIR, "disease", "donor_label")
    return library, donor, disease, taxonomy_lookup(cache)


def regional_metadata(cache, source_cells, chunk_size, library, donor, disease, tax_lookup):
    wanted = set(source_cells)
    cell_path = metadata_path(cache, DATASET_DIR, "cell_metadata")
    membership_path = metadata_path(cache, TAXONOMY_DIR, "cell_to_cluster_membership")

    logging.info("Filtering cell_metadata.csv for this region")
    cell = read_filtered_csv(cell_path, wanted, "cell_label", chunk_size)
    logging.info("Filtering cell_to_cluster_membership.csv for this region")
    membership = read_filtered_csv(
        membership_path, wanted, "cell_label", chunk_size,
        columns=["cell_label", "cluster_alias"],
    )
    if cell.empty or membership.empty:
        raise RuntimeError("No regional metadata rows were found")

    meta = cell.join(membership[["cluster_alias"]], how="inner")
    meta["cluster_alias"] = meta["cluster_alias"].astype("string")
    tax_lookup = tax_lookup.copy()
    tax_lookup.index = tax_lookup.index.astype("string")
    meta = meta.join(tax_lookup, on="cluster_alias", how="left")
    meta = meta.join(library, on="library_label", how="left", rsuffix="_library")
    donor_key = "donor_label_library" if "donor_label_library" in meta.columns else "donor_label"
    if donor_key not in meta.columns:
        raise KeyError("Regional metadata lacks donor_label")
    meta = meta.join(donor, on=donor_key, how="left", rsuffix="_donor")
    meta = meta.join(disease, on=donor_key, how="left", rsuffix="_disease")
    required = ["library_label", "cluster_alias", "Class", "Subclass", "Supertype"]
    missing = set(required).difference(meta.columns)
    if missing:
        raise KeyError(f"Regional metadata lacks {sorted(missing)}")
    meta = meta.dropna(subset=required)
    # Categories sharply reduce repeated-string storage within the region.
    for col in required + [c for c in ("Brain Region", "region_of_interest_label", "feature_matrix_label") if c in meta]:
        meta[col] = meta[col].astype("category")
    return meta


def matrix_names(cache, regions):
    available = cache.list_expression_matrix_files(directory=DATASET_DIR)
    available = {str(x).strip("/") for x in available}
    names = [f"{r}-10X/raw" for r in regions]
    missing = [x for x in names if x not in available]
    if missing:
        raise RuntimeError(f"Missing expression entries: {missing}")
    return names


def gene_metadata(cache, source):
    gene = read_small(cache, DATASET_DIR, "gene")
    gene["gene_identifier"] = gene["gene_identifier"].astype(str)
    gene = gene.drop_duplicates("gene_identifier").set_index("gene_identifier", drop=False)
    ids = pd.Index(source.var_names.astype(str), name="gene_identifier")
    var = gene.reindex(ids).copy()
    if var["gene_symbol"].isna().any():
        raise ValueError("Some source genes are absent from gene.csv")
    symbols = var["gene_symbol"].astype(str)
    counts = symbols.value_counts()
    keys = [s if counts[s] == 1 else f"{s}__{i}" for s, i in zip(symbols, ids)]
    var.insert(0, "gene_symbol_original", symbols.to_numpy())
    var.index = pd.Index(keys, name="gene_key")
    return var


def group_metadata(meta):
    frame = meta[["library_label", "Supertype"]].astype(str)
    mi = pd.MultiIndex.from_frame(frame)
    codes, groups = pd.factorize(mi, sort=True)
    keys = pd.DataFrame(groups.tolist(), columns=["library_label", "Supertype"])
    temp = meta.assign(_group_code=codes)
    first = temp.groupby("_group_code", observed=True, sort=True).first()
    first = first.drop(columns=["library_label", "Supertype"], errors="ignore")
    obs = pd.concat([keys, first.reset_index(drop=True)], axis=1)
    obs.insert(2, "number_of_nuclei", np.bincount(codes, minlength=len(obs)).astype(np.int64))
    obs.index = pd.Index(
        [f"{lib}__supertype_{alias}" for lib, alias in zip(obs["library_label"], obs["cluster_alias"].astype(str))],
        name="pseudobulk_id",
    )
    if not obs.index.is_unique:
        raise ValueError("Pseudobulk IDs are not unique")
    return obs, codes.astype(np.int64, copy=False)


def aggregate(source, full_codes, n_groups, chunk_size):
    counts = sparse.csr_matrix((n_groups, source.n_vars), dtype=np.int64)
    detected = sparse.csr_matrix((n_groups, source.n_vars), dtype=np.int64)
    for start in range(0, source.n_obs, chunk_size):
        stop = min(start + chunk_size, source.n_obs)
        local = full_codes[start:stop]
        keep = local >= 0
        if not keep.any():
            continue
        block = source.X[start:stop]
        block = block.tocsr() if sparse.issparse(block) else sparse.csr_matrix(np.asarray(block))
        block = block[keep]
        if block.nnz and not np.allclose(block.data, np.rint(block.data)):
            raise ValueError("Selected matrix is not integer raw counts")
        block.data = np.rint(block.data).astype(np.int64, copy=False)
        binary = block.copy()
        binary.data = np.ones(binary.nnz, dtype=np.int64)
        local = local[keep]
        assign = sparse.csr_matrix(
            (np.ones(len(local), dtype=np.int64), (local, np.arange(len(local)))),
            shape=(n_groups, len(local)),
        )
        counts = counts + assign @ block
        detected = detected + assign @ binary
        del block, binary, assign
        if start == 0 or stop == source.n_obs or (start // chunk_size) % 25 == 0:
            logging.info("Aggregated rows %d-%d of %d", start + 1, stop, source.n_obs)
    counts.eliminate_zeros(); detected.eliminate_zeros()
    return counts, detected.astype(np.uint32)


def sanitize(df):
    out = df.copy()
    for col in out.columns:
        if isinstance(out[col].dtype, pd.CategoricalDtype):
            out[col] = out[col].astype(str)
        elif out[col].dtype == object:
            out[col] = out[col].astype("string")
    return out


def app_obs_metadata(obs):
    """Return only the compact SEA-AD fields needed by the viewer.

    Besides fixing HDF5-incompatible source column names (for example names
    containing '/'), this prevents cell-level technical metadata from being
    repeated in every pseudobulk observation.
    """
    aliases = {
        "Supertype": ["Supertype"],
        "Brain Region": ["Brain Region", "region_of_interest_name", "region_of_interest_label", "feature_matrix_label"],
        "library_prep": ["library_prep", "library_method"],
        "Donor ID": ["Donor ID", "donor_label_library", "donor_label"],
        "Number of nuclei": ["Number of nuclei", "number_of_nuclei"],
        "Sex": ["Sex", "donor_sex", "donor_gender"],
        "Age at Death": ["Age at Death", "donor_age_value", "donor_age"],
        "Primary Study Name": ["Primary Study Name"],
        "Neurotypical reference": ["Neurotypical reference", "Neurotypical Reference"],
        "Overall AD neuropathological Change": ["Overall AD neuropathological Change"],
        "Thal": ["Thal"],
        "Braak": ["Braak"],
        "CERAD score": ["CERAD score"],
        "Cognitive Status": ["Cognitive Status"],
        "RIN": ["RIN"],
        "Severely Affected Donor": ["Severely Affected Donor"],
        "Brain pH": ["Brain pH"],
        "Class": ["Class"],
        "Subclass": ["Subclass"],
    }
    compact = pd.DataFrame(index=obs.index.copy())
    missing = []
    for target, candidates in aliases.items():
        source = next((name for name in candidates if name in obs.columns), None)
        if source is None:
            missing.append(target)
        else:
            compact[target] = obs[source].to_numpy()

    required = {"Supertype", "Brain Region", "Donor ID", "Number of nuclei", "Class", "Subclass"}
    missing_required = required.intersection(missing)
    if missing_required:
        raise KeyError(
            "Pseudobulk metadata lacks required viewer fields: "
            + ", ".join(sorted(missing_required))
        )
    if missing:
        logging.warning(
            "Optional viewer metadata fields not available and omitted: %s",
            ", ".join(missing),
        )
    compact["Number of nuclei"] = pd.to_numeric(
        compact["Number of nuclei"], errors="raise"
    ).astype(np.int64)
    return sanitize(compact)


def process_region(cache, matrix_name, scratch_dir, chunk_size, metadata_chunk_size,
                   library, donor, disease, tax_lookup, expected_var):
    source_path = Path(cache.get_data_path(directory=DATASET_DIR, file_name=matrix_name))
    logging.info("Opening %s in backed mode", source_path)
    source = ad.read_h5ad(source_path, backed="r")
    try:
        cells = pd.Index(source.obs_names.astype(str), name="cell_label")
        meta = regional_metadata(cache, cells, metadata_chunk_size, library, donor, disease, tax_lookup)
        aligned = meta.reindex(cells)
        assigned = aligned["library_label"].notna() & aligned["Supertype"].notna()
        logging.info("%s: %d/%d source nuclei assigned", matrix_name, int(assigned.sum()), source.n_obs)
        valid = aligned.loc[assigned].copy()
        obs, selected_codes = group_metadata(valid)
        full_codes = np.full(source.n_obs, -1, dtype=np.int64)
        full_codes[np.flatnonzero(assigned.to_numpy())] = selected_codes
        var = gene_metadata(cache, source)
        if expected_var is not None and not var.index.equals(expected_var.index):
            raise ValueError(f"Gene order differs in {matrix_name}")
        counts, detected = aggregate(source, full_codes, len(obs), chunk_size)
        compact_obs = app_obs_metadata(obs)
        region_file = scratch_dir / f"{matrix_name.split('/')[0]}_pseudobulk.h5ad"
        regional = ad.AnnData(
            X=counts, obs=compact_obs, var=sanitize(var),
            layers={"number_of_nuclei_expressing": detected},
        )
        regional.write_h5ad(region_file, compression="gzip")
        logging.info("Wrote regional checkpoint: %s", region_file)
        del meta, aligned, valid, obs, compact_obs, counts, detected, regional
        gc.collect()
        return region_file, var
    finally:
        source.file.close()
        del source
        gc.collect()


def main():
    a = args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    if not Path("/scratch").exists() and str(a.scratch_dir).startswith("/scratch"):
        raise RuntimeError("/scratch is not mounted; pass --scratch-dir to a writable large-volume path")
    a.cache_dir.mkdir(parents=True, exist_ok=True)
    a.scratch_dir.mkdir(parents=True, exist_ok=True)
    a.output.parent.mkdir(parents=True, exist_ok=True)

    cache = AbcProjectCache.from_cache_dir(a.cache_dir)
    load_manifest(cache, a.manifest)
    logging.info("Loading small lookup tables only")
    library, donor, disease, tax_lookup = static_metadata(cache)
    names = matrix_names(cache, a.regions)

    checkpoints, expected_var = [], None
    for name in names:
        expected_checkpoint = a.scratch_dir / f"{name.split('/')[0]}_pseudobulk.h5ad"
        if a.resume and expected_checkpoint.exists():
            logging.info("Reusing regional checkpoint: %s", expected_checkpoint)
            checkpoint_check = ad.read_h5ad(expected_checkpoint, backed="r")
            try:
                if "number_of_nuclei_expressing" not in checkpoint_check.layers:
                    raise RuntimeError(
                        f"Checkpoint lacks detection layer: {expected_checkpoint}"
                    )
                if "Number of nuclei" not in checkpoint_check.obs.columns:
                    raise RuntimeError(
                        f"Checkpoint lacks Number of nuclei: {expected_checkpoint}"
                    )
                checkpoint_var = checkpoint_check.var.copy()
            finally:
                checkpoint_check.file.close()
            if expected_var is not None and not checkpoint_var.index.equals(expected_var.index):
                raise ValueError(f"Gene order differs in checkpoint {expected_checkpoint}")
            checkpoint, var = expected_checkpoint, checkpoint_var
        else:
            checkpoint, var = process_region(
                cache, name, a.scratch_dir, a.matrix_chunk_size,
                a.metadata_chunk_size, library, donor, disease,
                tax_lookup, expected_var,
            )
        checkpoints.append(checkpoint)
        if expected_var is None:
            expected_var = var

    logging.info("Combining %d regional pseudobulk checkpoints", len(checkpoints))
    regional = [ad.read_h5ad(path) for path in checkpoints]
    if len(regional) == 1:
        combined = regional[0]
    else:
        # AnnData accepts only inner/outer joins. Gene identity and order were
        # already validated above, so an inner join is safe here.
        combined = ad.concat(
            regional, axis=0, join="inner", merge="same", index_unique=None
        )
    combined.uns.update({
        "dataset": DATASET_DIR,
        "taxonomy": TAXONOMY_DIR,
        "manifest": str(cache.current_manifest),
        "aggregation": "sum of raw counts by library_label and Supertype",
        "detection_layer": "number_of_nuclei_expressing",
        "fraction_expressing_formula": "number_of_nuclei_expressing / Number of nuclei",
        "regional_expression_matrices": names,
    })
    combined.write_h5ad(a.output, compression="gzip")
    check = ad.read_h5ad(a.output, backed="r")
    try:
        assert check.layers["number_of_nuclei_expressing"].shape == check.shape
        assert "Number of nuclei" in check.obs
    finally:
        check.file.close()
    if not a.keep_regional_files:
        for path in checkpoints:
            path.unlink(missing_ok=True)
    logging.info("Complete: %s rows x %s genes -> %s", combined.n_obs, combined.n_vars, a.output)


if __name__ == "__main__":
    main()
