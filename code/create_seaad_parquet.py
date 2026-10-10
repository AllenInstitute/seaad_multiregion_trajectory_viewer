#!/usr/bin/env python3
"""Convert the SEA-AD pseudobulk H5AD to one Parquet file per gene.

Supports resumable execution. With --resume, existing gene files are skipped only
when their schema and row count match the expected output contract.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
from pathlib import Path
from urllib.parse import quote

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse

DEFAULT_INPUT = Path(
    "/data/pseudobulk/SEA_AD_Multiregion_library_supertype_pseudobulk.h5ad"
)
DEFAULT_OUTPUT = Path("/results/SEAAD_parquet_by_gene")
DEFAULT_CPS_SOURCE = (
    "https://sea-ad-single-cell-profiling.s3.amazonaws.com/"
    "Multiregion_2026/previous_objects/parquet/metadata/"
    "global_cps_multiregion.csv"
)
DETECTION_LAYER = "number_of_nuclei_expressing"
NUCLEI_COLUMN = "Number of nuclei"
GENE_COLUMNS = (
    "sample_id",
    "summed_counts",
    "number_of_nuclei_expressing",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create one app-ready Parquet file per SEA-AD gene."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--gene-block-size", type=int, default=64)
    parser.add_argument(
        "--compression",
        default="zstd",
        choices=("zstd", "snappy", "gzip", "brotli", "none"),
    )
    parser.add_argument(
        "--cps-source",
        default=DEFAULT_CPS_SOURCE,
        help="Public URL or local CPS CSV path; use 'none' to skip.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete the existing output directory and start over.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Keep the existing output directory and skip valid completed "
            "gene Parquet files. Invalid or incomplete files are regenerated."
        ),
    )
    return parser.parse_args()


def prepare_output(path: Path, overwrite: bool, resume: bool) -> None:
    if overwrite and resume:
        raise ValueError("Use only one of --overwrite or --resume.")
    if path.exists():
        if overwrite:
            shutil.rmtree(path)
        elif not resume:
            raise FileExistsError(
                f"Output already exists: {path}. Use --resume to continue "
                "or --overwrite to replace it."
            )
    (path / "counts_by_gene").mkdir(parents=True, exist_ok=True)
    (path / "metadata").mkdir(parents=True, exist_ok=True)
    (path / "statistics").mkdir(parents=True, exist_ok=True)


def validate_input(adata: ad.AnnData) -> None:
    if DETECTION_LAYER not in adata.layers:
        raise KeyError(f"Missing layer: {DETECTION_LAYER}")
    if NUCLEI_COLUMN not in adata.obs.columns:
        raise KeyError(f"Missing obs column: {NUCLEI_COLUMN}")
    if adata.layers[DETECTION_LAYER].shape != adata.shape:
        raise ValueError("Detection layer shape does not match X.")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Input H5AD is empty.")
    if not adata.obs_names.is_unique:
        raise ValueError("Input observation names are not unique.")


def find_column(frame: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    return next((name for name in candidates if name in frame.columns), None)


def create_gene_map(var: pd.DataFrame) -> pd.DataFrame:
    symbol_column = find_column(
        var, ("gene_symbol", "gene_symbols", "symbol", "Gene Symbol")
    )
    identifiers = pd.Index(var.index.astype(str))
    if symbol_column is None:
        logging.warning("No gene-symbol column found; using var_names.")
        symbols = identifiers.copy()
    else:
        symbols = pd.Index(var[symbol_column].astype("string").fillna(""))
        if (symbols == "").any():
            values = symbols.to_numpy(dtype=object)
            missing = np.flatnonzero(values == "")
            values[missing] = identifiers.to_numpy()[missing]
            symbols = pd.Index(values.astype(str))

    seen: dict[str, int] = {}
    keys: list[str] = []
    for symbol in symbols.astype(str):
        base = quote(symbol, safe="") or "unnamed_gene"
        seen[base] = seen.get(base, 0) + 1
        keys.append(base if seen[base] == 1 else f"{base}__{seen[base]}")

    return pd.DataFrame(
        {
            "gene_symbol": symbols.astype(str),
            "gene_column": keys,
            "gene_identifier": identifiers.astype(str),
            "h5ad_column_index": np.arange(len(symbols), dtype=np.int32),
        }
    )


def prepare_obs(adata: ad.AnnData) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs.insert(0, "sample_id", adata.obs_names.astype(str))
    if obs["sample_id"].duplicated().any():
        raise ValueError("sample_id values are not unique.")

    # SEA-AD pseudobulk rows do not require a separate library_label field.
    # The app uses this column to verify exact row alignment with CPM_scaling.csv.
    # Therefore, when it is absent, use the stable unique pseudobulk sample_id.
    if "library_label" not in obs.columns:
        logging.info(
            "obs does not contain library_label; using sample_id as library_label."
        )
        obs["library_label"] = obs["sample_id"]
    else:
        obs["library_label"] = obs["library_label"].astype("string")
        missing = obs["library_label"].isna() | (
            obs["library_label"].str.strip() == ""
        )
        if missing.any():
            logging.warning(
                "Replacing %d missing library_label values with sample_id.",
                int(missing.sum()),
            )
            obs.loc[missing, "library_label"] = obs.loc[missing, "sample_id"]
    obs[NUCLEI_COLUMN] = pd.to_numeric(obs[NUCLEI_COLUMN], errors="raise")
    if obs[NUCLEI_COLUMN].isna().any() or (obs[NUCLEI_COLUMN] <= 0).any():
        raise ValueError("Number of nuclei must be finite and positive.")
    for column in obs.columns:
        if isinstance(obs[column].dtype, pd.CategoricalDtype):
            obs[column] = obs[column].astype("string")
    return obs.reset_index(drop=True)


def create_cpm_scaling(adata: ad.AnnData, obs: pd.DataFrame) -> pd.DataFrame:
    """Match the mouse pipeline: CPM uses each pseudobulk row's total raw counts."""
    # Slicing a backed CSRDataset loads a scipy sparse matrix supporting sum.
    row_totals = np.zeros(adata.n_obs, dtype=np.float64)
    row_block_size = 1024
    matrix = adata.X
    for start in range(0, adata.n_obs, row_block_size):
        stop = min(start + row_block_size, adata.n_obs)
        block = matrix[start:stop, :]
        if sparse.issparse(block):
            totals = np.asarray(block.sum(axis=1, dtype=np.float64)).ravel()
        else:
            totals = np.asarray(block, dtype=np.float64).sum(axis=1)
        row_totals[start:stop] = totals
        del block
    logging.info("Calculated CPM library totals for %d pseudobulk rows.", adata.n_obs)
    if row_totals.shape[0] != adata.n_obs:
        raise RuntimeError("Unexpected row count while calculating CPM scaling.")
    if not np.all(np.isfinite(row_totals)) or np.any(row_totals < 0):
        raise ValueError("Pseudobulk library totals must be finite and nonnegative.")

    scaling = np.zeros_like(row_totals, dtype=np.float64)
    positive = row_totals > 0
    scaling[positive] = 1_000_000.0 / row_totals[positive]

    return pd.DataFrame(
        {
            "library_label": obs["library_label"].astype(str).to_numpy(),
            "CPM_scaling_factor": scaling,
        }
    )


def matrix_block(matrix, start: int, stop: int):
    block = matrix[:, start:stop]
    return block.tocsc() if sparse.issparse(block) else np.asarray(block)


def extract_column(block, offset: int) -> np.ndarray:
    if sparse.issparse(block):
        return block.getcol(offset).toarray().ravel()
    return np.asarray(block[:, offset]).ravel()


def integer_counts(block, offset: int, dtype) -> np.ndarray:
    values = extract_column(block, offset)
    if not np.all(np.isfinite(values)):
        raise ValueError("Matrix contains nonfinite values.")
    if np.any(values < 0) or not np.allclose(values, np.rint(values)):
        raise ValueError("Matrix contains negative or noninteger values.")
    rounded = np.rint(values)
    if np.issubdtype(dtype, np.unsignedinteger):
        if rounded.size and rounded.max() > np.iinfo(dtype).max:
            raise OverflowError(f"Values exceed {dtype} capacity.")
    return rounded.astype(dtype, copy=False)


def valid_gene_file(path: Path, expected_rows: int) -> bool:
    """Validate metadata only, avoiding a full Parquet read during resume."""
    if not path.is_file():
        return False
    try:
        parquet_file = pq.ParquetFile(path)
        return (
            parquet_file.metadata.num_rows == expected_rows
            and set(parquet_file.schema_arrow.names) == set(GENE_COLUMNS)
        )
    except Exception as error:
        logging.warning("Existing file is invalid and will be replaced: %s (%s)", path, error)
        return False


def write_gene_files(
    adata: ad.AnnData,
    gene_map: pd.DataFrame,
    output: Path,
    block_size: int,
    compression: str,
    resume: bool,
) -> None:
    sample_ids = np.asarray(adata.obs_names.astype(str), dtype=object)
    nuclei = pd.to_numeric(adata.obs[NUCLEI_COLUMN], errors="raise").to_numpy(
        dtype=np.int64
    )
    parquet_compression = None if compression == "none" else compression
    written = 0
    skipped = 0

    for start in range(0, adata.n_vars, block_size):
        stop = min(start + block_size, adata.n_vars)
        pending: list[int] = []

        for gene_index in range(start, stop):
            key = gene_map.iloc[gene_index]["gene_column"]
            path = output / "counts_by_gene" / f"gene={key}" / "data.parquet"
            if resume and valid_gene_file(path, adata.n_obs):
                skipped += 1
            else:
                pending.append(gene_index)

        if not pending:
            logging.info("Skipping completed genes %d-%d of %d", start + 1, stop, adata.n_vars)
            continue

        logging.info(
            "Reading genes %d-%d of %d (%d need writing)",
            start + 1,
            stop,
            adata.n_vars,
            len(pending),
        )
        counts_block = matrix_block(adata.X, start, stop)
        detection_block = matrix_block(adata.layers[DETECTION_LAYER], start, stop)

        for gene_index in pending:
            offset = gene_index - start
            counts = integer_counts(counts_block, offset, np.int64)
            detected = integer_counts(detection_block, offset, np.uint32)
            if np.any(detected.astype(np.int64) > nuclei):
                symbol = gene_map.iloc[gene_index]["gene_symbol"]
                raise ValueError(
                    f"Detection counts exceed Number of nuclei for gene {symbol!r}."
                )

            key = gene_map.iloc[gene_index]["gene_column"]
            directory = output / "counts_by_gene" / f"gene={key}"
            directory.mkdir(parents=True, exist_ok=True)
            destination = directory / "data.parquet"
            temporary = directory / "data.parquet.tmp"
            temporary.unlink(missing_ok=True)

            table = pa.table(
                {
                    "sample_id": pa.array(sample_ids, type=pa.string()),
                    "summed_counts": pa.array(counts, type=pa.int64()),
                    "number_of_nuclei_expressing": pa.array(
                        detected, type=pa.uint32()
                    ),
                }
            )
            pq.write_table(
                table,
                temporary,
                compression=parquet_compression,
                use_dictionary=["sample_id"],
                write_statistics=True,
            )
            temporary.replace(destination)
            written += 1

        del counts_block, detection_block

    logging.info(
        "Gene-file summary: wrote %d; skipped %d valid existing files",
        written,
        skipped,
    )


def preserve_cps(source: str, metadata_dir: Path) -> None:
    if source.strip().lower() == "none":
        logging.info("Skipping CPS metadata copy")
        return
    logging.info("Reading CPS metadata from %s", source)
    cps = pd.read_csv(source, low_memory=False)
    required = [
        "donor_id",
        "CPS_Global",
        "CPS_Global_ABeta",
        "CPS_Global_pTau",
    ]
    missing = [column for column in required if column not in cps.columns]
    if missing:
        raise KeyError(f"CPS file is missing columns: {missing}")
    cps = cps[required].copy()
    if cps["donor_id"].isna().any() or cps["donor_id"].duplicated().any():
        raise ValueError("CPS donor_id values must be present and unique.")
    for column in required[1:]:
        cps[column] = pd.to_numeric(cps[column], errors="raise")
    cps.to_csv(metadata_dir / "global_cps_multiregion.csv", index=False)


def validate_output(
    output: Path,
    gene_map: pd.DataFrame,
    n_obs: int,
    expected_library_labels: np.ndarray,
) -> None:
    if pq.ParquetFile(output / "metadata" / "obs_metadata.parquet").metadata.num_rows != n_obs:
        raise RuntimeError("obs_metadata.parquet row-count mismatch.")
    cpm_path = output / "metadata" / "CPM_scaling.csv"
    saved_cpm = pd.read_csv(cpm_path)
    if list(saved_cpm.columns) != ["library_label", "CPM_scaling_factor"]:
        raise RuntimeError("CPM_scaling.csv columns do not match the app contract.")
    if len(saved_cpm) != n_obs:
        raise RuntimeError("CPM_scaling.csv row-count mismatch.")
    if not np.array_equal(
        saved_cpm["library_label"].astype(str).to_numpy(),
        expected_library_labels.astype(str),
    ):
        raise RuntimeError("CPM_scaling.csv library order does not match metadata.")
    scaling = pd.to_numeric(saved_cpm["CPM_scaling_factor"], errors="coerce")
    if scaling.isna().any() or (scaling < 0).any():
        raise RuntimeError("CPM_scaling.csv contains invalid scaling factors.")

    saved_map = pq.read_table(output / "metadata" / "gene_map.parquet").to_pandas()
    if len(saved_map) != len(gene_map):
        raise RuntimeError("gene_map.parquet row-count mismatch.")

    invalid = []
    for key in saved_map["gene_column"]:
        path = output / "counts_by_gene" / f"gene={key}" / "data.parquet"
        if not valid_gene_file(path, n_obs):
            invalid.append(str(key))
            if len(invalid) >= 10:
                break
    if invalid:
        raise RuntimeError(f"Missing or invalid gene files, including: {invalid}")


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    if args.gene_block_size < 1:
        raise ValueError("--gene-block-size must be at least 1.")
    if not args.input.is_file():
        raise FileNotFoundError(f"Input H5AD not found: {args.input}")

    prepare_output(args.output, args.overwrite, args.resume)
    logging.info("Opening pseudobulk H5AD in backed mode: %s", args.input)
    adata = ad.read_h5ad(args.input, backed="r")
    try:
        validate_input(adata)
        obs = prepare_obs(adata)
        gene_map = create_gene_map(adata.var.copy())
        cpm_scaling = create_cpm_scaling(adata, obs)
        compression = None if args.compression == "none" else args.compression

        # These small files are safe to rewrite on a resumed run.
        pq.write_table(
            pa.Table.from_pandas(obs, preserve_index=False),
            args.output / "metadata" / "obs_metadata.parquet",
            compression=compression,
        )
        pq.write_table(
            pa.Table.from_pandas(gene_map, preserve_index=False),
            args.output / "metadata" / "gene_map.parquet",
            compression=compression,
        )
        cpm_scaling.to_csv(
            args.output / "metadata" / "CPM_scaling.csv",
            index=False,
        )

        write_gene_files(
            adata,
            gene_map,
            args.output,
            args.gene_block_size,
            args.compression,
            args.resume,
        )
        preserve_cps(args.cps_source, args.output / "metadata")

        manifest = {
            "input_h5ad": str(args.input),
            "n_pseudobulk_samples": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "gene_file_layout": "counts_by_gene/gene=<gene_column>/data.parquet",
            "gene_columns": list(GENE_COLUMNS),
            "cpm_definition": (
                "summed_counts * (1,000,000 / total raw counts across all genes "
                "in the pseudobulk row)"
            ),
            "sample_order": (
                "Identical across obs_metadata.parquet, CPM_scaling.csv, and "
                "every gene file"
            ),
            "fraction_expressing": (
                "number_of_nuclei_expressing / Number of nuclei"
            ),
        }
        with open(args.output / "manifest.json", "w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2)

        validate_output(
            args.output,
            gene_map,
            adata.n_obs,
            obs["library_label"].astype(str).to_numpy(),
        )
        logging.info("Completed SEA-AD Parquet conversion: %s", args.output)
    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
