#!/usr/bin/env python3
"""Plan a SEA-AD partition-key migration without recomputing expression.
Requires anndata, pandas, numpy, pyarrow. AWS execution uses the AWS CLI.
Run against the ORIGINAL pseudobulk H5AD and its corresponding original gene map.
Outputs new metadata, an audit CSV, and an executable AWS migration script.
Does not modify input files or execute AWS operations.
"""
import argparse
import shlex
from pathlib import Path
from urllib.parse import quote
import anndata as ad
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--h5ad', required=True, type=Path)
    p.add_argument('--gene-map', required=True, type=Path)
    p.add_argument('--s3-root', required=True, help='S3 prefix containing counts_by_gene and metadata')
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--statistics', nargs='*', type=Path, default=[],
                   help='Optional local CSV/Parquet statistics files with gene_column keys; published under statistics/')
    a = p.parse_args()
    if not a.s3_root.startswith('s3://'):
        raise ValueError('--s3-root must be an S3 URI')
    if a.output.exists():
        raise FileExistsError('Use a new output directory for each migration plan.')
    old = pd.read_parquet(a.gene_map)
    required = {'gene_symbol', 'gene_column', 'h5ad_column_index'}
    if not required.issubset(old.columns):
        raise ValueError(f'Gene map must contain {sorted(required)}')
    h = ad.read_h5ad(a.h5ad, backed='r')
    try:
        var = h.var.copy()
    finally:
        h.file.close()
    if len(old) != len(var):
        raise ValueError('H5AD and gene-map gene counts differ.')
    idx_values = pd.to_numeric(old['h5ad_column_index'], errors='raise').to_numpy()
    if not np.all(np.isfinite(idx_values)) or not np.all(idx_values == np.floor(idx_values)):
        raise ValueError('Invalid H5AD column indices.')
    idx = idx_values.astype(int)
    if set(idx) != set(range(len(var))):
        raise ValueError('H5AD column indices must be a permutation of all gene columns.')
    v = var.iloc[idx]
    for c in ('gene_symbol', 'gene_identifier'):
        if c not in v or v[c].isna().any() or (v[c].astype(str).str.strip() == '').any():
            raise ValueError(f'Original H5AD must have complete var[{c!r}].')
    symbols = v['gene_symbol'].astype(str).to_numpy()
    if not np.array_equal(old['gene_symbol'].astype(str).to_numpy(), symbols):
        raise ValueError('Gene symbols do not align: wrong H5AD/map pair.')
    # Reconstruct original SEA-AD Parquet keys in source-column order.
    seen = {}
    original_keys = []
    for symbol in var['gene_symbol'].astype(str):
        base = quote(symbol, safe='') or 'unnamed_gene'
        seen[base] = seen.get(base, 0) + 1
        original_keys.append(base if seen[base] == 1 else f'{base}__{seen[base]}')
    oldkeys = old['gene_column'].astype(str).to_numpy()
    if not np.array_equal(oldkeys, np.asarray(original_keys)[idx]):
        raise ValueError('Map does not match original SEA-AD naming; do not migrate an already changed map.')
    newkeys = np.asarray([quote(str(k), safe='') for k in v.index])
    if len(set(oldkeys)) != len(oldkeys) or len(set(newkeys)) != len(newkeys):
        raise ValueError('Old or new partition keys collide. No plan generated.')
    duplicate = pd.Series(symbols).duplicated(keep=False).to_numpy()
    changed = oldkeys != newkeys
    if np.any(changed & ~duplicate):
        raise ValueError('A unique-symbol key would change; inspect the H5AD var_names.')
    new = old.copy()
    new['gene_column'] = newkeys
    new['gene_identifier'] = v['gene_identifier'].astype(str).to_numpy()
    if 'gene_key' in new:
        new['gene_key'] = newkeys
    mapping = dict(zip(oldkeys, newkeys))
    # Prepare optional statistics updates before writing anything.
    updates = []
    names = set()
    for path in a.statistics:
        if path.name in names:
            raise ValueError('Statistics basenames must be unique.')
        names.add(path.name)
        if path.suffix.lower() == '.csv':
            frame = pd.read_csv(path)
        elif path.suffix.lower() == '.parquet':
            frame = pd.read_parquet(path)
        else:
            raise ValueError('Statistics must be CSV or Parquet.')
        if 'gene_column' not in frame:
            raise ValueError(f'{path}: no gene_column; inspect its schema separately.')
        if frame['gene_column'].isna().any() or not frame['gene_column'].astype(str).isin(mapping).all():
            raise ValueError(f'{path}: unknown or missing gene keys.')
        frame['gene_column'] = frame['gene_column'].astype(str).map(mapping)
        if 'gene_key' in frame:
            if not frame['gene_key'].astype(str).isin(mapping).all():
                raise ValueError(f'{path}: unknown gene_key values.')
            frame['gene_key'] = frame['gene_key'].astype(str).map(mapping)
        updates.append((path.name, frame))
    a.output.mkdir(parents=True)
    newpath = (a.output / 'gene_map.parquet').resolve()
    new.to_parquet(newpath, index=False)
    audit = pd.DataFrame({'gene_symbol': symbols, 'old_key': oldkeys, 'new_key': newkeys,
                          'changed': changed, 'gene_identifier': new['gene_identifier']})
    audit.to_csv(a.output / 'gene_key_migration.csv', index=False)
    for name, frame in updates:
        target = a.output / name
        if target.suffix.lower() == '.csv':
            frame.to_csv(target, index=False)
        else:
            frame.to_parquet(target, index=False)
    root = a.s3_root.rstrip('/')
    q = shlex.quote
    lines = ['#!/usr/bin/env bash', 'set -euo pipefail',
             '# Run only after reviewing the audit CSV. AWS credentials must be configured.',
             '# All changed sources are staged first to avoid overwrite/name-overlap corruption.',
             'STAGE=' + q(root + '/_gene_key_migration_') + '"$(date -u +%Y%m%dT%H%M%SZ)_${RANDOM}"',
             'echo "Staging under $STAGE"']
    pairs = [(str(o), str(n)) for o,n,c in zip(oldkeys,newkeys,changed) if c]
    for i,(o,n) in enumerate(pairs):
        source = root + '/counts_by_gene/gene=' + o + '/data.parquet'
        lines.append(f'aws s3 cp {q(source)} "$STAGE/{i}/data.parquet" --only-show-errors')
    for i,(o,n) in enumerate(pairs):
        dest = root + '/counts_by_gene/gene=' + n + '/data.parquet'
        lines.append(f'aws s3 cp "$STAGE/{i}/data.parquet" {q(dest)} --only-show-errors')
    # AWS CLI copy success is required before publishing the replacement map.
    for name,frame in updates:
        lines.append(f'aws s3 cp {q(str((a.output/name).resolve()))} {q(root+"/statistics/"+name)} --only-show-errors')
    lines.append(f'aws s3 cp {q(str(newpath))} {q(root+"/metadata/gene_map.parquet")} --only-show-errors')
    # Never delete a key that is also a destination in the final map.
    finalkeys = set(newkeys)
    for o,n in pairs:
        if o not in finalkeys:
            lines.append(f'aws s3 rm {q(root+"/counts_by_gene/gene="+o+"/data.parquet")} --only-show-errors')
    lines += ['aws s3 rm "$STAGE/" --recursive --only-show-errors',
              'echo "Migration complete. Restart the app to clear cached gene maps."']
    shell = a.output / 'apply_aws_migration.sh'
    shell.write_text('\n'.join(lines)+'\n')
    shell.chmod(0o755)
    print(f'Genes: {len(old)}; renamed partitions: {int(changed.sum())}; unique-symbol partitions unchanged.')
    print(f'Output: {a.output.resolve()}')
    print('Review gene_key_migration.csv, then run bash apply_aws_migration.sh.')
    print('Only explicitly supplied statistics files are updated; inspect other key-dependent resources.')

if __name__ == '__main__':
    main()
