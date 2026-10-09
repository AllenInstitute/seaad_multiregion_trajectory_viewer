#!/usr/bin/env python3
"""Create SEA-AD categorical value sets from pseudobulk observation metadata.

The output contains the columns required by the Shiny app:
``label,field,color_hex_triplet,order``.

Taxonomy colors/orders are taken from the pseudobulk H5AD when available.
Clinically meaningful stage orders are defined explicitly. Remaining observed
categorical values receive deterministic colors and first-observed order.
Continuous fields, including CPS values, are intentionally excluded.
"""

from __future__ import annotations

import argparse
import colorsys
import hashlib
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

DEFAULT_INPUT = "data/pseudobulk/SEA_AD_Multiregion_library_supertype_pseudobulk.h5ad"
DEFAULT_OUTPUT = "data/parquet/metadata/value_sets.csv"

ALIASES = {
    "Brain Region": ["Brain Region", "region_of_interest_label", "feature_matrix_label"],
    "Sex": ["Sex", "donor_sex", "donor_gender"],
    "library_prep": ["library_prep", "library_method"],
}

EXPLICIT_ORDERS = {
    "Overall AD neuropathological Change": ["Not AD", "Low", "Intermediate", "High"],
    "Thal": ["Thal 0", "Thal 1", "Thal 2", "Thal 3", "Thal 4", "Thal 5"],
    "Braak": ["Braak 0", "Braak I", "Braak II", "Braak III", "Braak IV", "Braak V", "Braak VI"],
    "CERAD score": ["Absent", "Sparse", "Moderate", "Frequent"],
    "Cognitive Status": ["No dementia", "Dementia"],
    "Severely Affected Donor": ["N", "Y"],
}

CANDIDATE_FIELDS = [
    "Class", "Subclass", "Supertype", "Brain Region", "library_prep",
    "Sex", "Primary Study Name", "Neurotypical reference",
    "Neurotypical Reference", "Overall AD neuropathological Change", "Thal",
    "Braak", "CERAD score", "Cognitive Status", "Severely Affected Donor",
    "APOE Genotype", "Overall CAA Score", "Highest Lewy Body Disease", "LATE",
    "Atherosclerosis", "Arteriolosclerosis",
]

CONTINUOUS_FIELDS = {
    "CPS_Global", "CPS_Global_ABeta", "CPS_Global_pTau", "Age at Death",
    "Brain pH", "RIN", "Number of nuclei", "number_of_nuclei",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=Path(DEFAULT_INPUT))
    p.add_argument("--output", type=Path, default=Path(DEFAULT_OUTPUT))
    return p.parse_args()


def deterministic_color(field: str, label: str) -> str:
    digest = hashlib.sha1(f"{field}\0{label}".encode()).digest()
    hue = int.from_bytes(digest[:2], "big") / 65535.0
    saturation = 0.45 + digest[2] / 255.0 * 0.25
    lightness = 0.42 + digest[3] / 255.0 * 0.20
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return f"#{round(r*255):02X}{round(g*255):02X}{round(b*255):02X}"


def first_present(columns, candidates):
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def normalized_metadata(obs: pd.DataFrame) -> pd.DataFrame:
    result = obs.copy()
    for target, candidates in ALIASES.items():
        source = first_present(result.columns, candidates)
        if source is not None and target not in result.columns:
            result[target] = result[source]
    return result


def ordered_values(obs: pd.DataFrame, field: str) -> list[str]:
    observed = obs[field].dropna().astype(str)
    observed = [value for value in pd.unique(observed) if value and value.lower() != "nan"]
    specified = EXPLICIT_ORDERS.get(field, [])
    return [x for x in specified if x in observed] + [x for x in observed if x not in specified]


def taxonomy_spec(obs: pd.DataFrame, field: str):
    color_column = f"{field}_color"
    order_column = f"{field}_order"
    if color_column not in obs.columns and order_column not in obs.columns:
        return {}
    columns = [field]
    if color_column in obs.columns:
        columns.append(color_column)
    if order_column in obs.columns:
        columns.append(order_column)
    frame = obs[columns].dropna(subset=[field]).copy()
    frame[field] = frame[field].astype(str)
    if order_column in frame.columns:
        frame[order_column] = pd.to_numeric(frame[order_column], errors="coerce")
    frame = frame.drop_duplicates(field)
    result = {}
    for _, row in frame.iterrows():
        result[row[field]] = {
            "color": row.get(color_column),
            "order": row.get(order_column),
        }
    return result


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    if not args.input.exists():
        raise FileNotFoundError(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input, backed="r")
    try:
        obs = normalized_metadata(adata.obs)
    finally:
        adata.file.close()

    rows = []
    for field in CANDIDATE_FIELDS:
        if field in CONTINUOUS_FIELDS or field not in obs.columns:
            continue
        values = ordered_values(obs, field)
        taxonomy = taxonomy_spec(obs, field)
        ranked = []
        for fallback_order, label in enumerate(values, start=1):
            saved_order = taxonomy.get(label, {}).get("order")
            order = int(saved_order) if pd.notna(saved_order) else fallback_order
            color = taxonomy.get(label, {}).get("color")
            if not isinstance(color, str) or not color.startswith("#") or len(color) != 7:
                color = deterministic_color(field, label)
            ranked.append({
                "label": label,
                "field": field,
                "color_hex_triplet": color,
                "order": order,
            })
        rows.extend(sorted(ranked, key=lambda row: (row["order"], row["label"])))

    result = pd.DataFrame(rows, columns=["label", "field", "color_hex_triplet", "order"])
    if result.empty:
        raise RuntimeError("No categorical value sets were generated")
    if result.duplicated(["field", "label"]).any():
        raise RuntimeError("Duplicate field/label rows were generated")
    result.to_csv(args.output, index=False)

    check = pd.read_csv(args.output)
    required = {"label", "field", "color_hex_triplet", "order"}
    if not required.issubset(check.columns):
        raise RuntimeError("Written value_sets.csv lacks required columns")
    logging.info("Complete: %s", args.output)


if __name__ == "__main__":
    main()
