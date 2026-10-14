#!/usr/bin/env python3
"""Merge reference metadata and taxonomy value sets without replacing source colors/orders.

Required reference inputs are explicit to avoid accidentally reading the generated output.
Optional --input H5AD adds observed categorical labels absent from the references.
All reference labels are retained, including labels not observed in the H5AD.
Taxonomy is authoritative for Class, Subclass and Supertype.
Only genuinely missing observed labels receive deterministic fallback colors/orders.
"""
from __future__ import annotations
import argparse
import colorsys
import csv
import hashlib
import logging
import io
from urllib.request import urlopen
from urllib.parse import urlparse
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path

COLUMNS = ['label', 'field', 'color_hex_triplet', 'order']
ALIASES = {
    'region_of_interest_label': 'Brain Region',
    'region_of_interest': 'Brain Region',
    'brain_region': 'Brain Region',
    'donor_sex': 'Sex',
    'donor_gender': 'Sex',
    'library_method': 'library_prep',
    'class': 'Class', 'subclass': 'Subclass', 'supertype': 'Supertype',
}
CANDIDATE_FIELDS = [
    'Class', 'Subclass', 'Supertype', 'Brain Region', 'library_prep', 'Sex',
    'Primary Study Name', 'Neurotypical reference', 'Neurotypical Reference',
    'Overall AD neuropathological Change', 'Thal', 'Braak', 'CERAD score',
    'Cognitive Status', 'Severely Affected Donor', 'APOE Genotype',
    'Overall CAA Score', 'Highest Lewy Body Disease', 'LATE',
    'Atherosclerosis', 'Arteriolosclerosis', 'dementia',
]

def canonical_field(value):
    value = str(value).strip()
    return ALIASES.get(value.lower(), value)

def is_url(source):
    return urlparse(str(source)).scheme in {'http', 'https'}

def read_csv(path):
    if is_url(path):
        with urlopen(str(path), timeout=120) as response:
            text = response.read().decode('utf-8-sig')
        handle = io.StringIO(text, newline='')
    else:
        handle = Path(path).open(encoding='utf-8-sig', newline='')
    with handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        rows = list(reader)
    return headers, rows

def validated(label, field, color, order, context):
    label, field = str(label).strip(), canonical_field(field)
    color, order = str(color).strip(), str(order).strip()
    if not label or not field:
        raise ValueError(f'{context}: empty label or field')
    if not re.fullmatch(r'#[0-9a-fA-F]{6}', color):
        raise ValueError(f'{context}: invalid color {color!r}; source not silently replaced')
    try:
        number = Decimal(order)
        if not number.is_finite():
            raise ValueError('non-finite order')
    except (InvalidOperation, ValueError) as exc:
        raise ValueError(f'{context}: invalid order {order!r}') from exc
    # Preserve source spelling, case and order values, including gaps.
    return dict(zip(COLUMNS, [label, field, color, order]))

def add_unique(target, row, context):
    key = (row['field'], row['label'])
    previous = target.get(key)
    if previous is not None:
        if (previous['color_hex_triplet'].lower() != row['color_hex_triplet'].lower()
                or Decimal(previous['order']) != Decimal(row['order'])):
            raise ValueError(f'{context}: conflicting definitions for {key}')
    else:
        target[key] = row

def deterministic_color(field, label):
    digest = hashlib.sha1(f'{field}\0{label}'.encode()).digest()
    hue = int.from_bytes(digest[:2], 'big') / 65535.0
    saturation = 0.45 + digest[2] / 255.0 * 0.25
    lightness = 0.42 + digest[3] / 255.0 * 0.20
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return f'#{round(r*255):02X}{round(g*255):02X}{round(b*255):02X}'

def build(metadata_path, taxonomy_path, h5ad_path=None):
    headers, source = read_csv(metadata_path)
    missing = set(COLUMNS) - set(headers)
    if missing:
        raise ValueError(f'{metadata_path}: missing columns {sorted(missing)}')
    metadata = {}
    for index, row in enumerate(source, 2):
        add_unique(metadata, validated(*(row[c] for c in COLUMNS),
                   f'{metadata_path}:{index}'), str(metadata_path))
    headers, source = read_csv(taxonomy_path)
    required = {'name', 'cluster_annotation_term_set_name', 'color_hex_triplet', 'term_order'}
    if required - set(headers):
        raise ValueError(f'{taxonomy_path}: missing columns {sorted(required-set(headers))}. Use cluster_annotation_term.csv, not cluster_annotation_term_set.csv.')
    taxonomy = {}
    for index, row in enumerate(source, 2):
        field = canonical_field(row['cluster_annotation_term_set_name'])
        if field not in {'Class', 'Subclass', 'Supertype'}:
            continue
        # Display names, NOT accession labels; term_order, NOT row position or term_set_order.
        item = validated(row['name'], field, row['color_hex_triplet'], row['term_order'],
                         f'{taxonomy_path}:{index}')
        add_unique(taxonomy, item, str(taxonomy_path))
    if not taxonomy:
        raise ValueError('No Class/Subclass/Supertype taxonomy rows found')
    result = dict(metadata)
    for key, row in taxonomy.items():
        if key in result and result[key] != row:
            logging.warning('Taxonomy overrides metadata definition for %s', key)
        result[key] = row
    reference = dict(result)
    if h5ad_path is not None:
        import anndata as ad
        import pandas as pd
        adata = ad.read_h5ad(h5ad_path, backed='r')
        try:
            obs = adata.obs.copy()
        finally:
            adata.file.close()
        # Gather aliases together; never discard one column's values in favor of another.
        observed = {}
        relevant = set(CANDIDATE_FIELDS) | {field for field, _ in result}
        for column in obs.columns:
            field = canonical_field(column)
            if field not in relevant:
                continue
            labels = {str(v).strip() for v in obs[column] if not pd.isna(v) and str(v).strip()}
            observed.setdefault(field, set()).update(labels)
        for field, labels in observed.items():
            missing_labels = sorted(label for label in labels if (field, label) not in result)
            orders = [Decimal(row['order']) for row in result.values() if row['field'] == field]
            next_order = int(max(orders, default=Decimal(0)).to_integral_value(rounding='ROUND_FLOOR')) + 1
            for label in missing_labels:
                logging.warning('No reference definition for %s / %s; assigning fallback', field, label)
                result[(field, label)] = validated(label, field, deterministic_color(field, label),
                                                  str(next_order), 'fallback')
                next_order += 1
    for key, row in reference.items():
        if result[key] != row:
            raise AssertionError(f'Reference definition was modified: {key}')
    return sorted(result.values(), key=lambda row: (row['field'], Decimal(row['order']), row['label']))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--value-sets', '--metadata-value-sets', dest='value_sets', required=True)
    parser.add_argument('--taxonomy', '--cluster-annotation-term', dest='taxonomy', required=True)
    parser.add_argument('--input', type=Path, default=None, help='Optional pseudobulk H5AD for missing observed labels')
    parser.add_argument('--output', type=Path, default=Path('data/parquet/metadata/value_sets.csv'))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s | %(message)s')
    for source in [args.value_sets, args.taxonomy, args.input]:
        if source is not None and not is_url(source):
            source = Path(source)
            if not source.is_file():
                raise FileNotFoundError(source)
            if source.resolve() == args.output.resolve():
                raise ValueError('Output must not overwrite an input reference file')
    rows = build(args.value_sets, args.taxonomy, args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + '.tmp')
    try:
        with temporary.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        headers, check = read_csv(temporary)
        if headers != COLUMNS or check != rows:
            raise RuntimeError('Written output did not match validated rows')
        temporary.replace(args.output)
    finally:
        temporary.unlink(missing_ok=True)
    logging.info('Wrote %s rows to %s', len(rows), args.output)

if __name__ == '__main__':
    main()
