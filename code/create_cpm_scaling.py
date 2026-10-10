#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


input_h5ad = Path(
    "../data/pseudobulk/"
    "SEA_AD_Multiregion_library_supertype_pseudobulk.h5ad"
)

output_csv = Path(
    "../scrattch/SEAAD_parquet_by_gene/metadata/CPM_scaling.csv"
)


# Load the complete H5AD into memory. This avoids issues with summing
# backed sparse datasets.
adata = ad.read_h5ad(input_h5ad)

try:
    # Each observation must uniquely represent one library-supertype row.
    if not adata.obs_names.is_unique:
        duplicated_names = (
            pd.Index(adata.obs_names)[
                pd.Index(adata.obs_names).duplicated(keep=False)
            ]
            .unique()
            .tolist()
        )

        raise ValueError(
            "H5AD obs_names must be unique. "
            f"First duplicated names: {duplicated_names[:10]}"
        )

    row_label = adata.obs_names.astype(str).to_numpy()

    # Extract the library identifier from observation names such as:
    # L8XR_231012_01_G04__supertype_175
    separator = "__supertype_"

    valid_name_format = np.char.find(
        row_label.astype(str),
        separator,
    ) > 0

    if not np.all(valid_name_format):
        bad_rows = np.flatnonzero(~valid_name_format)

        bad_names = row_label[bad_rows[:10]].tolist()

        raise ValueError(
            "One or more obs_names do not follow the expected format "
            "'<library_label>__supertype_<supertype_id>'. "
            f"First affected names: {bad_names}"
        )

    library_label = np.asarray(
        [
            label.split(separator, 1)[0]
            for label in row_label
        ],
        dtype=str,
    )

    supertype_id = np.asarray(
        [
            label.split(separator, 1)[1]
            for label in row_label
        ],
        dtype=str,
    )

    # Store the reconstructed library label in obs for consistency.
    adata.obs["library_label"] = library_label

    if np.any(library_label == ""):
        bad_rows = np.flatnonzero(library_label == "")

        raise ValueError(
            "One or more reconstructed library labels are empty. "
            f"First affected rows: {bad_rows[:10].tolist()}"
        )

    if np.any(supertype_id == ""):
        bad_rows = np.flatnonzero(supertype_id == "")

        raise ValueError(
            "One or more reconstructed supertype IDs are empty. "
            f"First affected rows: {bad_rows[:10].tolist()}"
        )

    # Calculate the total count for each library-supertype pseudobulk row.
    # Summing before conversion preserves sparse-matrix behavior.
    if sparse.issparse(adata.X):
        total_counts = np.asarray(
            adata.X.sum(axis=1)
        ).ravel()
    else:
        total_counts = np.asarray(
            adata.X.sum(axis=1)
        ).ravel()

    total_counts = np.asarray(
        total_counts,
        dtype=np.float64,
    )

    if total_counts.shape[0] != adata.n_obs:
        raise ValueError(
            "The number of calculated row totals does not match n_obs. "
            f"Calculated {total_counts.shape[0]} totals for "
            f"{adata.n_obs} observations."
        )

    if not np.all(np.isfinite(total_counts)):
        bad_rows = np.flatnonzero(
            ~np.isfinite(total_counts)
        )

        raise ValueError(
            "The pseudobulk matrix produced non-finite row totals. "
            f"First affected rows: {bad_rows[:10].tolist()}"
        )

    if np.any(total_counts <= 0):
        bad_rows = np.flatnonzero(total_counts <= 0)

        bad_names = row_label[bad_rows[:10]].tolist()

        raise ValueError(
            "One or more pseudobulk observations have zero or negative "
            "total counts. "
            f"First affected rows: {bad_rows[:10].tolist()}; "
            f"observation names: {bad_names}"
        )

    cpm_scaling = pd.DataFrame(
        {
            "row_label": row_label,
            "library_label": library_label,
            "supertype_id": supertype_id,
            "total_counts": total_counts,
            "CPM_scaling_factor": 1_000_000.0 / total_counts,
        }
    )

    # library_label is intentionally non-unique because each library has
    # multiple supertype pseudobulk rows. row_label must remain unique.
    if not cpm_scaling["row_label"].is_unique:
        raise ValueError(
            "row_label must be unique for row-aligned app metadata."
        )

    if len(cpm_scaling) != adata.n_obs:
        raise ValueError(
            "The CPM scaling table is not aligned with the H5AD rows."
        )

    output_csv.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    cpm_scaling.to_csv(
        output_csv,
        index=False,
    )

finally:
    # Close the AnnData file manager if one is present.
    if getattr(adata, "file", None) is not None:
        adata.file.close()


print(f"Wrote {output_csv}")
print(f"Rows written: {len(cpm_scaling):,}")
print(
    f"Unique libraries: "
    f"{cpm_scaling['library_label'].nunique():,}"
)
print()
print(cpm_scaling.head())
print()
print(cpm_scaling["CPM_scaling_factor"].describe())
