from pathlib import Path
import sys
import pyarrow.parquet as pq

root = Path("data/parquet")
counts_dir = root / "counts_by_gene"
obs_file = root / "metadata" / "obs_metadata.parquet"
gene_map_file = root / "metadata" / "gene_map.parquet"

required = [counts_dir, obs_file, gene_map_file]
missing = [str(path) for path in required if not path.exists()]

if missing:
    sys.exit(f"Missing required paths: {missing}")

# Validate metadata files.
obs = pq.read_table(obs_file)
gene_map = pq.read_table(gene_map_file).to_pandas()

if "sample_id" not in obs.column_names:
    sys.exit(
        "obs_metadata.parquet lacks sample_id. "
        f"Columns: {obs.column_names}"
    )

sample_ids = obs["sample_id"].to_pylist()

if len(sample_ids) != len(set(sample_ids)):
    sys.exit("sample_id is not unique in obs_metadata.parquet")

# Locate the gene-partition key in gene_map.parquet.
possible_key_columns = [
    "gene_column",
    "gene_key",
    "partition",
    "gene_partition",
]

key_col = next(
    (c for c in possible_key_columns if c in gene_map.columns),
    None,
)

if key_col is None:
    sys.exit(
        "Could not identify the gene partition key in gene_map.parquet.\n"
        f"Columns found: {gene_map.columns.tolist()}"
    )

expected_keys = set(gene_map[key_col].dropna().astype(str))

partition_dirs = sorted(
    path for path in counts_dir.iterdir()
    if path.is_dir()
)

def partition_key(path):
    return path.name.split("=", 1)[-1]

observed_keys = {partition_key(path) for path in partition_dirs}

missing_partitions = sorted(expected_keys - observed_keys)
unexpected_partitions = sorted(observed_keys - expected_keys)

empty_partitions = []
invalid_partitions = []
total_files = 0
expected_rows = len(sample_ids)

for index, directory in enumerate(partition_dirs, start=1):
    files = sorted(directory.rglob("*.parquet"))
    total_files += len(files)

    if not files:
        empty_partitions.append(str(directory))
        continue

    rows = 0
    columns = set()

    try:
        for file in files:
            parquet = pq.ParquetFile(file)
            rows += parquet.metadata.num_rows
            columns.update(parquet.schema_arrow.names)
    except Exception as exc:
        invalid_partitions.append(
            (str(directory), f"unreadable: {exc}")
        )
        continue

    required_columns = {"sample_id", "summed_counts"}

    if not required_columns.issubset(columns):
        invalid_partitions.append(
            (str(directory), f"columns={sorted(columns)}")
        )
    elif rows != expected_rows:
        invalid_partitions.append(
            (
                str(directory),
                f"{rows} rows; expected {expected_rows}",
            )
        )

    if index % 1000 == 0:
        print(
            f"Checked {index:,}/{len(partition_dirs):,} partitions"
        )

print()
print(f"Observation rows:       {expected_rows:,}")
print(f"Rows in gene map:       {len(gene_map):,}")
print(f"Expected gene keys:     {len(expected_keys):,}")
print(f"Gene directories:       {len(partition_dirs):,}")
print(f"Parquet fragments:      {total_files:,}")
print(f"Missing partitions:     {len(missing_partitions):,}")
print(f"Unexpected partitions:  {len(unexpected_partitions):,}")
print(f"Empty partitions:       {len(empty_partitions):,}")
print(f"Invalid partitions:     {len(invalid_partitions):,}")

if missing_partitions:
    print("First missing:", missing_partitions[:20])

if unexpected_partitions:
    print("First unexpected:", unexpected_partitions[:20])

if empty_partitions:
    print("First empty:", empty_partitions[:20])

if invalid_partitions:
    print("First invalid:", invalid_partitions[:20])

if (
    missing_partitions
    or unexpected_partitions
    or empty_partitions
    or invalid_partitions
):
    sys.exit("PARQUET VALIDATION FAILED")

print("FULL PARQUET VALIDATION PASSED")
