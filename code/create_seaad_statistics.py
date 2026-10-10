#!/usr/bin/env python3
"""Create the SEA-AD gene summary table using the mouse-app calculations.

The calculations match the developing-mouse implementation, except that the
age comparison is replaced by Dementia versus No dementia. Expression values
are CPM calculated from each pseudobulk library total, not counts per nucleus.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from pybiomart import Dataset

DEFAULT_INPUT = Path(
    "data/pseudobulk/SEA_AD_Multiregion_library_supertype_pseudobulk.h5ad"
)
DEFAULT_OUTPUT = Path(
    "/results/SEAAD_parquet_by_gene/statistics/seaad_gene_summary_statistics.csv"
)
GENE_BLOCK_SIZE = 256
MIN_NONZERO_N = 3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--gene-block-size", type=int, default=GENE_BLOCK_SIZE)
    return parser.parse_args()


def first_existing(columns, candidates):
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def infer_gene_symbols(adata):
    for field in (
        "gene_symbol_original",
        "gene_symbol",
        "gene_symbols",
        "symbol",
        "gene_name",
    ):
        if field in adata.var.columns:
            values = adata.var[field].astype("string")
            fallback = pd.Series(
                adata.var_names.astype(str),
                index=adata.var_names,
                dtype="string",
            )
            return (
                values.replace("", pd.NA)
                .fillna(fallback)
                .astype(str)
                .to_numpy()
            )
    return adata.var_names.astype(str).to_numpy()


def fetch_gene_types(gene_symbols):
    import io
    import requests

    gene_symbols = np.asarray(
        [str(value).strip() for value in gene_symbols],
        dtype=object,
    )

    query_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" '
        'header="1" uniqueRows="1" datasetConfigVersion="0.6">'
        '<Dataset name="hsapiens_gene_ensembl" interface="default">'
        '<Attribute name="external_gene_name"/>'
        '<Attribute name="gene_biotype"/>'
        '</Dataset>'
        '</Query>'
    )

    hosts = [
        "https://www.ensembl.org",
        "https://useast.ensembl.org",
        "https://asia.ensembl.org",
    ]

    df = None
    for host in hosts:
        url = f"{host}/biomart/martservice"
        try:
            logging.info("Querying BioMart: %s", url)
            response = requests.get(
                url,
                params={"query": query_xml},
                timeout=300,
            )
            response.raise_for_status()
            text = response.text

            if text.lstrip().lower().startswith("<"):
                logging.warning("BioMart at %s returned HTML; trying next host.", host)
                continue

            candidate = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
            candidate.columns = ["gene_symbol", "gene_biotype"]

            if candidate.empty:
                logging.warning("BioMart at %s returned no rows; trying next host.", host)
                continue

            df = candidate
            logging.info("Retrieved %d BioMart rows from %s", len(df), host)
            break

        except Exception as error:
            logging.warning("BioMart query failed at %s: %s", host, error)

    if df is None:
        logging.warning(
            "All BioMart hosts failed; gene_type will be 'unknown' for all genes."
        )
        return np.full(len(gene_symbols), "unknown", dtype=object)

    df = df.dropna(subset=["gene_symbol", "gene_biotype"])
    df["gene_symbol"] = df["gene_symbol"].str.strip()
    df = df[df["gene_symbol"] != ""].drop_duplicates(subset=["gene_symbol"])

    gene_type_map = dict(zip(df["gene_symbol"], df["gene_biotype"]))

    unmatched = [s for s in dict.fromkeys(gene_symbols) if s not in gene_type_map]
    logging.info(
        "BioMart matched %d of %d unique SEA-AD gene symbols",
        len(set(gene_symbols)) - len(unmatched),
        len(set(gene_symbols)),
    )
    if unmatched:
        logging.info("Example unmatched gene symbols: %s", unmatched[:20])

    return np.asarray(
        [gene_type_map.get(symbol, "unknown") for symbol in gene_symbols],
        dtype=object,
    )

def normalize_dementia_status(values):
    """Map Cognitive Status to the app's binary dementia field."""
    text = pd.Series(values, dtype="string").str.strip()
    lower = text.str.lower()
    result = pd.Series(pd.NA, index=text.index, dtype="string")

    no_dementia = lower.str.contains(
        r"no dementia|not demented|cognitively normal",
        regex=True,
        na=False,
    )
    result.loc[no_dementia] = "No dementia"
    result.loc[text.notna() & (text != "") & ~no_dementia] = "Dementia"
    return result.to_numpy(dtype=object)


def safe_log2fc(numerator, denominator):
    return np.log2(
        (np.asarray(numerator, dtype=np.float64) + 1.0)
        / (np.asarray(denominator, dtype=np.float64) + 1.0)
    )


def positive_mean(values, minimum_n=MIN_NONZERO_N):
    positive = values > 0
    positive_n = positive.sum(axis=0)
    positive_sum = np.where(positive, values, 0.0).sum(axis=0)
    result = np.full(values.shape[1], np.nan, dtype=np.float64)
    eligible = positive_n >= minimum_n
    result[eligible] = positive_sum[eligible] / positive_n[eligible]
    return result


def grouped_means(
    values,
    groups,
    nonzero_only=False,
    minimum_nonzero_n=MIN_NONZERO_N,
):
    groups = np.asarray(groups, dtype=object)
    valid_group = np.asarray(
        [
            value is not None
            and not pd.isna(value)
            and str(value).strip() != ""
            for value in groups
        ],
        dtype=bool,
    )
    group_strings = np.asarray(
        ["" if not valid else str(value) for value, valid in zip(groups, valid_group)],
        dtype=object,
    )
    labels = pd.unique(group_strings[valid_group]).tolist()
    means = []
    for label in labels:
        group_values = values[valid_group & (group_strings == label), :]
        if nonzero_only:
            group_mean = positive_mean(
                group_values,
                minimum_n=minimum_nonzero_n,
            )
        else:
            group_mean = np.mean(group_values, axis=0)
        means.append(group_mean)
    if not means:
        return [], np.empty((0, values.shape[1]), dtype=np.float64)
    return labels, np.vstack(means)


def named_group_mean(labels, means, label):
    if label not in labels:
        return np.full(means.shape[1], np.nan, dtype=np.float64)
    return means[labels.index(label), :]


def maximum_group_statistics(labels, means):
    n_genes = means.shape[1]
    maximum_labels = np.empty(n_genes, dtype=object)
    maximum_labels[:] = None
    log2fc = np.full(n_genes, np.nan, dtype=np.float64)
    for gene_index in range(n_genes):
        gene_means = means[:, gene_index]
        valid = np.isfinite(gene_means)
        if not np.any(valid):
            continue
        valid_indices = np.flatnonzero(valid)
        maximum_index = valid_indices[np.argmax(gene_means[valid])]
        maximum_value = gene_means[maximum_index]
        mean_across_groups = np.mean(gene_means[valid])
        maximum_labels[gene_index] = labels[maximum_index]
        log2fc[gene_index] = safe_log2fc(maximum_value, mean_across_groups)
    return maximum_labels, log2fc


def weighted_fraction(detected_values, cell_totals):
    cell_totals = np.asarray(cell_totals, dtype=np.float64)
    denominator = np.sum(cell_totals)
    if not np.isfinite(denominator) or denominator <= 0:
        return np.full(detected_values.shape[1], np.nan, dtype=np.float64)
    return np.sum(detected_values, axis=0, dtype=np.float64) / denominator


def grouped_detection_fractions(detected_values, cell_totals, groups):
    groups = np.asarray(groups, dtype=object)
    cell_totals = np.asarray(cell_totals, dtype=np.float64)
    valid_group = np.asarray(
        [
            value is not None
            and not pd.isna(value)
            and str(value).strip() != ""
            for value in groups
        ],
        dtype=bool,
    )
    group_strings = np.asarray(
        ["" if not valid else str(value) for value, valid in zip(groups, valid_group)],
        dtype=object,
    )
    labels = pd.unique(group_strings[valid_group]).tolist()
    fractions = []
    for label in labels:
        mask = valid_group & (group_strings == label)
        fractions.append(
            weighted_fraction(detected_values[mask, :], cell_totals[mask])
        )
    if not fractions:
        return [], np.empty((0, detected_values.shape[1]), dtype=np.float64)
    return labels, np.vstack(fractions)


def maximum_fraction_statistics(labels, fractions):
    n_genes = fractions.shape[1]
    maximum_labels = np.empty(n_genes, dtype=object)
    maximum_labels[:] = None
    fraction_difference = np.full(n_genes, np.nan, dtype=np.float64)
    for gene_index in range(n_genes):
        gene_fractions = fractions[:, gene_index]
        valid = np.isfinite(gene_fractions)
        if not np.any(valid):
            continue
        valid_indices = np.flatnonzero(valid)
        maximum_index = valid_indices[np.argmax(gene_fractions[valid])]
        maximum_labels[gene_index] = labels[maximum_index]
        fraction_difference[gene_index] = (
            gene_fractions[maximum_index] - np.mean(gene_fractions[valid])
        )
    return maximum_labels, fraction_difference


def dense_block(matrix, start, stop):
    block = matrix[:, start:stop]
    if sparse.issparse(block):
        return block.toarray().astype(np.float64, copy=False)
    return np.asarray(block, dtype=np.float64)


def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    if not args.input.exists():
        candidates = list(Path("/data").glob("**/*pseudobulk*.h5ad"))
        if len(candidates) == 1:
            args.input = candidates[0]
        else:
            raise FileNotFoundError(args.input)

    logging.info("Reading: %s", args.input)
    adata = ad.read_h5ad(args.input)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"Input H5AD is empty: {adata.shape}")

    cognitive_field = first_existing(
        adata.obs.columns,
        ("Cognitive Status", "dementia", "cognitive_status"),
    )
    roi_field = first_existing(
        adata.obs.columns,
        (
            "Brain Region",
            "region_of_interest_label",
            "brain_region",
            "region_label",
            "structure",
            "ROI",
        ),
    )
    subclass_field = first_existing(
        adata.obs.columns,
        ("Subclass", "subclass", "subclass_label"),
    )
    cell_count_field = first_existing(
        adata.obs.columns,
        ("Number of nuclei", "number_of_cells", "number_of_nuclei"),
    )

    if cognitive_field is None:
        raise KeyError("Could not identify Cognitive Status/dementia column.")
    if roi_field is None:
        raise KeyError("Could not identify the brain-region column.")
    if subclass_field is None:
        raise KeyError("Could not identify the subclass column.")
    if cell_count_field is None:
        raise KeyError("Could not identify the number-of-cells/nuclei column.")

    if cognitive_field == "dementia":
        dementia_groups = adata.obs[cognitive_field].astype("string").to_numpy(dtype=object)
    else:
        dementia_groups = normalize_dementia_status(adata.obs[cognitive_field])

    observed_statuses = set(pd.Series(dementia_groups).dropna().astype(str))
    required_statuses = {"Dementia", "No dementia"}
    if not required_statuses.issubset(observed_statuses):
        raise ValueError(
            "Dementia comparison requires both 'Dementia' and 'No dementia'; "
            f"found {sorted(observed_statuses)}"
        )

    roi_groups = adata.obs[roi_field].astype("string").to_numpy(dtype=object)
    subclass_groups = adata.obs[subclass_field].astype("string").to_numpy(dtype=object)
    number_of_cells = pd.to_numeric(
        adata.obs[cell_count_field], errors="coerce"
    ).to_numpy(dtype=np.float64)
    if np.any(~np.isfinite(number_of_cells)) or np.any(number_of_cells <= 0):
        raise ValueError("Number of cells/nuclei must contain finite positive values.")

    gene_symbols = infer_gene_symbols(adata)
    gene_types = fetch_gene_types(gene_symbols)

    expression_matrix = adata.X
    detection_layer = first_existing(
        list(adata.layers.keys()),
        ("number_of_cells_expressing", "number_of_nuclei_expressing"),
    )
    if detection_layer is None:
        raise KeyError(
            "Input H5AD must contain number_of_cells_expressing or "
            "number_of_nuclei_expressing."
        )
    detection_matrix = adata.layers[detection_layer]
    if detection_matrix.shape != adata.shape:
        raise ValueError("Detection-count layer shape does not match X.")

    library_totals = np.asarray(expression_matrix.sum(axis=1)).reshape(-1).astype(np.float64)
    if np.any(~np.isfinite(library_totals)) or np.any(library_totals <= 0):
        raise ValueError("Every pseudobulk library must have a finite positive count total.")
    cpm_scaling_factor = 1_000_000.0 / library_totals

    statistics_blocks = []
    for block_start in range(0, adata.n_vars, args.gene_block_size):
        block_stop = min(block_start + args.gene_block_size, adata.n_vars)
        raw_counts = dense_block(expression_matrix, block_start, block_stop)
        detected_counts = dense_block(detection_matrix, block_start, block_stop)

        if np.any(detected_counts < 0):
            raise ValueError("Detection counts contain negative values.")
        if np.any(detected_counts > number_of_cells[:, np.newaxis]):
            raise ValueError(
                "Detection counts exceed the number of cells/nuclei for at least one row."
            )

        cpm = raw_counts * cpm_scaling_factor[:, np.newaxis]
        ln_mean_cpm = np.log1p(np.mean(cpm, axis=0))

        dementia_labels, dementia_means = grouped_means(
            cpm, dementia_groups, nonzero_only=False
        )
        dementia_mean = named_group_mean(dementia_labels, dementia_means, "Dementia")
        no_dementia_mean = named_group_mean(
            dementia_labels, dementia_means, "No dementia"
        )
        dementia_log2fc = safe_log2fc(dementia_mean, no_dementia_mean)

        dementia_nz_labels, dementia_nz_means = grouped_means(
            cpm, dementia_groups, nonzero_only=True
        )
        dementia_nz_mean = named_group_mean(
            dementia_nz_labels, dementia_nz_means, "Dementia"
        )
        no_dementia_nz_mean = named_group_mean(
            dementia_nz_labels, dementia_nz_means, "No dementia"
        )
        dementia_log2fc_nonzero = safe_log2fc(
            dementia_nz_mean, no_dementia_nz_mean
        )

        roi_labels, roi_means = grouped_means(cpm, roi_groups)
        max_roi, roi_log2fc = maximum_group_statistics(roi_labels, roi_means)
        roi_nz_labels, roi_nz_means = grouped_means(
            cpm, roi_groups, nonzero_only=True
        )
        _, roi_log2fc_nonzero = maximum_group_statistics(
            roi_nz_labels, roi_nz_means
        )

        subclass_labels, subclass_means = grouped_means(cpm, subclass_groups)
        max_subclass, subclass_log2fc = maximum_group_statistics(
            subclass_labels, subclass_means
        )
        subclass_nz_labels, subclass_nz_means = grouped_means(
            cpm, subclass_groups, nonzero_only=True
        )
        _, subclass_log2fc_nonzero = maximum_group_statistics(
            subclass_nz_labels, subclass_nz_means
        )

        overall_fraction = weighted_fraction(detected_counts, number_of_cells)
        dementia_fraction_labels, dementia_fractions = grouped_detection_fractions(
            detected_counts, number_of_cells, dementia_groups
        )
        dementia_fraction = named_group_mean(
            dementia_fraction_labels, dementia_fractions, "Dementia"
        )
        no_dementia_fraction = named_group_mean(
            dementia_fraction_labels, dementia_fractions, "No dementia"
        )
        dementia_fraction_difference = dementia_fraction - no_dementia_fraction

        roi_fraction_labels, roi_fractions = grouped_detection_fractions(
            detected_counts, number_of_cells, roi_groups
        )
        max_fraction_roi, fraction_difference_roi = maximum_fraction_statistics(
            roi_fraction_labels, roi_fractions
        )
        subclass_fraction_labels, subclass_fractions = grouped_detection_fractions(
            detected_counts, number_of_cells, subclass_groups
        )
        max_fraction_subclass, fraction_difference_subclass = (
            maximum_fraction_statistics(
                subclass_fraction_labels, subclass_fractions
            )
        )

        statistics_blocks.append(
            pd.DataFrame(
                {
                    "gene_symbol": gene_symbols[block_start:block_stop],
                    "ln(mean_CPM+1)": ln_mean_cpm,
                    "dementia_log2fc": dementia_log2fc,
                    "dementia_log2fc(>0)": dementia_log2fc_nonzero,
                    "max_ROI": max_roi,
                    "log2fc_ROI": roi_log2fc,
                    "log2fc_ROI(>0)": roi_log2fc_nonzero,
                    "max_subclass": max_subclass,
                    "log2fc_subclass": subclass_log2fc,
                    "log2fc_subclass(>0)": subclass_log2fc_nonzero,
                    "gene_type": gene_types[block_start:block_stop],
                    "frac_expr": overall_fraction,
                    "dementia_frac_diff": dementia_fraction_difference,
                    "max_frac_ROI": max_fraction_roi,
                    "frac_diff_ROI": fraction_difference_roi,
                    "max_frac_subclass": max_fraction_subclass,
                    "frac_diff_subclass": fraction_difference_subclass,
                }
            )
        )
        logging.info(
            "Completed genes %d-%d of %d",
            block_start + 1,
            block_stop,
            adata.n_vars,
        )

    statistics = pd.concat(statistics_blocks, ignore_index=True)
    expected_columns = [
        "gene_symbol",
        "ln(mean_CPM+1)",
        "dementia_log2fc",
        "dementia_log2fc(>0)",
        "max_ROI",
        "log2fc_ROI",
        "log2fc_ROI(>0)",
        "max_subclass",
        "log2fc_subclass",
        "log2fc_subclass(>0)",
        "gene_type",
        "frac_expr",
        "dementia_frac_diff",
        "max_frac_ROI",
        "frac_diff_ROI",
        "max_frac_subclass",
        "frac_diff_subclass",
    ]
    statistics = statistics[expected_columns]

    if statistics.shape[0] != adata.n_vars:
        raise RuntimeError("Statistics table does not contain one row per gene.")
    if statistics["gene_symbol"].isna().any():
        raise RuntimeError("Statistics table contains missing gene symbols.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    statistics.to_csv(args.output, index=False)

    check = pd.read_csv(args.output)
    if list(check.columns) != expected_columns:
        raise RuntimeError("Output column names changed during serialization.")
    if check.shape[0] != adata.n_vars:
        raise RuntimeError("Output row count changed during serialization.")

    logging.info("Saved: %s", args.output)
    logging.info("Rows: %d; columns: %d", *statistics.shape)


if __name__ == "__main__":
    main()
