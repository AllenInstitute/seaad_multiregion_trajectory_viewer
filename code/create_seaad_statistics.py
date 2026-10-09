#!/usr/bin/env python3
"""Create SEA-AD gene summary statistics from the pseudobulk H5AD.

The established SEA-AD table contains disease, brain-region and cell-type
summaries. This implementation calculates each summary twice:

* using all pseudobulk observations;
* using only pseudobulks with nonzero expression for the gene.

Expression is counts per 10,000 nuclei, matching the Shiny app:
``summed_counts / Number of nuclei * 10000``.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

DEFAULT_INPUT = "data/pseudobulk/SEA_AD_Multiregion_library_supertype_pseudobulk.h5ad"
DEFAULT_OUTPUT = "/results/SEAAD_parquet_by_gene/statistics/SEAAD_gene_summary_statistics.csv"
EPSILON = 1e-8


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=Path(DEFAULT_INPUT))
    p.add_argument("--output", type=Path, default=Path(DEFAULT_OUTPUT))
    p.add_argument("--gene-block-size", type=int, default=64)
    p.add_argument("--skip-fraction-statistics", action="store_true",
                   help="Do not calculate statistics from expressing-nucleus fractions.")
    return p.parse_args()


def first_present(columns, values):
    for value in values:
        if value in columns:
            return value
    return None


def gene_symbols(adata: ad.AnnData) -> np.ndarray:
    if "gene_symbol_original" in adata.var.columns:
        return adata.var["gene_symbol_original"].astype(str).to_numpy()
    if "gene_symbol" in adata.var.columns:
        return adata.var["gene_symbol"].astype(str).to_numpy()
    return adata.var_names.astype(str).to_numpy()


def matrix_block(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def mean_for_mask(values: np.ndarray, mask: np.ndarray) -> float:
    selected = values[mask & np.isfinite(values)]
    return float(np.mean(selected)) if selected.size else np.nan


def log2_fc(a: float, b: float) -> float:
    return float(np.log2((a + EPSILON) / (b + EPSILON))) if np.isfinite(a) and np.isfinite(b) else np.nan


def two_group_summary(values, groups, positive_only=False):
    valid = np.isfinite(values) & pd.notna(groups)
    if positive_only:
        valid &= values > 0
    g = np.asarray(groups, dtype=object)
    dementia = mean_for_mask(values, valid & (g == "Dementia"))
    no_dementia = mean_for_mask(values, valid & (g == "No dementia"))
    statistic = np.nan
    p_value = np.nan
    x = values[valid & (g == "Dementia")]
    y = values[valid & (g == "No dementia")]
    if len(x) >= 2 and len(y) >= 2:
        test = stats.mannwhitneyu(x, y, alternative="two-sided")
        statistic, p_value = float(test.statistic), float(test.pvalue)
    return dementia, no_dementia, log2_fc(dementia, no_dementia), statistic, p_value


def specificity_summary(values, labels, positive_only=False):
    labels = np.asarray(labels, dtype=object)
    valid = np.isfinite(values) & pd.notna(labels)
    if positive_only:
        valid &= values > 0
    if not np.any(valid):
        return None, np.nan, np.nan, np.nan
    frame = pd.DataFrame({"value": values[valid], "label": labels[valid]})
    means = frame.groupby("label", observed=True)["value"].mean().sort_values(ascending=False)
    if means.empty:
        return None, np.nan, np.nan, np.nan
    winner = str(means.index[0])
    maximum = float(means.iloc[0])
    other = frame.loc[frame["label"] != means.index[0], "value"]
    others_mean = float(other.mean()) if len(other) else np.nan
    return winner, maximum, others_mean, log2_fc(maximum, others_mean)


def spearman_summary(values, progression):
    progression = pd.to_numeric(pd.Series(progression), errors="coerce").to_numpy(float)
    valid = np.isfinite(values) & np.isfinite(progression)
    if valid.sum() < 3 or np.nanstd(values[valid]) == 0 or np.nanstd(progression[valid]) == 0:
        return np.nan, np.nan, int(valid.sum())
    result = stats.spearmanr(progression[valid], values[valid])
    return float(result.statistic), float(result.pvalue), int(valid.sum())


def summarize_gene(symbol, expression, obs, fraction=None):
    cognitive = obs["Cognitive Status"].astype(str).to_numpy() if "Cognitive Status" in obs else np.full(len(obs), None)
    region = obs["Brain Region"].astype(str).to_numpy()
    subclass = obs["Subclass"].astype(str).to_numpy()
    row = {"gene_symbol": symbol}

    for suffix, positive in (("", False), ("_nonzero", True)):
        d, n, lfc, u, p = two_group_summary(expression, cognitive, positive)
        row.update({
            f"mean_expression_Dementia{suffix}": d,
            f"mean_expression_No_dementia{suffix}": n,
            f"dementia_log2_fold_change{suffix}": lfc,
            f"dementia_mann_whitney_u{suffix}": u,
            f"dementia_p_value{suffix}": p,
        })
        winner, top, rest, spec = specificity_summary(expression, region, positive)
        row.update({
            f"maximum_expression_Brain_Region{suffix}": winner,
            f"maximum_Brain_Region_mean{suffix}": top,
            f"other_Brain_Region_mean{suffix}": rest,
            f"Brain_Region_specificity_log2_fold_change{suffix}": spec,
        })
        winner, top, rest, spec = specificity_summary(expression, subclass, positive)
        row.update({
            f"maximum_expression_Subclass{suffix}": winner,
            f"maximum_Subclass_mean{suffix}": top,
            f"other_Subclass_mean{suffix}": rest,
            f"Subclass_specificity_log2_fold_change{suffix}": spec,
        })

    for progression_name in ("CPS_Global", "CPS_Global_ABeta", "CPS_Global_pTau"):
        if progression_name in obs.columns:
            rho, p, n = spearman_summary(expression, obs[progression_name])
            row[f"{progression_name}_Spearman_rho"] = rho
            row[f"{progression_name}_Spearman_p_value"] = p
            row[f"{progression_name}_n"] = n

    if fraction is not None:
        for suffix, positive in (("", False), ("_nonzero", True)):
            d, n, lfc, u, p = two_group_summary(fraction, cognitive, positive)
            row.update({
                f"mean_fraction_expressing_Dementia{suffix}": d,
                f"mean_fraction_expressing_No_dementia{suffix}": n,
                f"dementia_fraction_expressing_log2_fold_change{suffix}": lfc,
                f"dementia_fraction_mann_whitney_u{suffix}": u,
                f"dementia_fraction_p_value{suffix}": p,
            })
            winner, top, rest, spec = specificity_summary(fraction, region, positive)
            row.update({
                f"maximum_fraction_expressing_Brain_Region{suffix}": winner,
                f"Brain_Region_fraction_specificity_log2_fold_change{suffix}": spec,
            })
            winner, top, rest, spec = specificity_summary(fraction, subclass, positive)
            row.update({
                f"maximum_fraction_expressing_Subclass{suffix}": winner,
                f"Subclass_fraction_specificity_log2_fold_change{suffix}": spec,
            })
    return row


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    if not args.input.exists():
        raise FileNotFoundError(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input, backed="r")
    try:
        obs = adata.obs.copy()
        count_column = first_present(obs.columns, ["number_of_nuclei", "Number of nuclei", "number_of_cells"])
        if count_column is None:
            raise KeyError("Input H5AD lacks number_of_nuclei/Number of nuclei")
        region_column = first_present(obs.columns, ["Brain Region", "region_of_interest_label", "feature_matrix_label"])
        subclass_column = first_present(obs.columns, ["Subclass", "subclass"])
        if region_column is None or subclass_column is None:
            raise KeyError("Input H5AD lacks Brain Region or Subclass")
        obs = obs.rename(columns={region_column: "Brain Region", subclass_column: "Subclass"})
        nuclei = pd.to_numeric(obs[count_column], errors="coerce").to_numpy(float)
        if np.any(~np.isfinite(nuclei)) or np.any(nuclei <= 0):
            raise ValueError("Nucleus counts must be finite and positive")

        symbols = gene_symbols(adata)
        detection_layer = first_present(
            list(adata.layers.keys()),
            ["number_of_nuclei_expressing", "number_of_cells_expressing"],
        )
        rows = []
        for start in range(0, adata.n_vars, args.gene_block_size):
            stop = min(start + args.gene_block_size, adata.n_vars)
            logging.info("Summarizing genes %d-%d of %d", start + 1, stop, adata.n_vars)
            counts = matrix_block(adata.X[:, start:stop]).astype(float, copy=False)
            expression = counts / nuclei[:, None] * 10000.0
            fractions = None
            if detection_layer is not None and not args.skip_fraction_statistics:
                detected = matrix_block(adata.layers[detection_layer][:, start:stop]).astype(float, copy=False)
                fractions = detected / nuclei[:, None]
            for local, gene_index in enumerate(range(start, stop)):
                rows.append(summarize_gene(
                    symbols[gene_index],
                    expression[:, local],
                    obs,
                    None if fractions is None else fractions[:, local],
                ))
        result = pd.DataFrame(rows)
        result.to_csv(args.output, index=False)
    finally:
        adata.file.close()

    check = pd.read_csv(args.output)
    if "gene_symbol" not in check.columns or len(check) == 0:
        raise RuntimeError("Written statistics table is invalid")
    logging.info("Complete: %s", args.output)


if __name__ == "__main__":
    main()
