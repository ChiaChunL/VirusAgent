#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

# Normalize columns and parse bait/strain from filenames

def parse_variant_bait(path: Path):
    name = path.stem  # e.g., Alpha-NSP12 or Alpha-CON
    if '-' in name:
        variant, bait = name.split('-', 1)
    else:
        variant, bait = 'Unknown', name
    is_control = bait.upper() == 'CON'
    return variant, bait, is_control


def read_apms_table(xlsx: Path) -> pd.DataFrame:
    df = pd.read_excel(xlsx)
    cols = {c: c for c in df.columns}
    # Unify name column
    if 'ProteinName' in cols:
        df.rename(columns={'ProteinName': 'Protein names'}, inplace=True)
    # Find intensity column (starts with 'Intensity')
    intensity_col = next((c for c in df.columns if str(c).startswith('Intensity')), None)
    if intensity_col is None:
        raise ValueError(f'No Intensity column in {xlsx}')
    df = df[['Protein IDs', 'Protein names', 'Fasta headers', intensity_col]].copy()
    df.rename(columns={intensity_col: 'intensity'}, inplace=True)
    # Add metadata
    variant, bait, is_control = parse_variant_bait(xlsx)
    df['strain'] = variant
    df['bait'] = bait
    df['is_control'] = is_control
    df['source_file'] = str(xlsx)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apms_root', required=True)
    ap.add_argument('--out_csv', required=True)
    args = ap.parse_args()

    root = Path(args.apms_root)
    rows = []
    for sub in root.glob('*'):
        if sub.is_dir():
            for xlsx in sub.glob('*.xlsx'):
                rows.append(read_apms_table(xlsx))
        elif sub.suffix.lower() == '.xlsx':
            rows.append(read_apms_table(sub))
    long_df = pd.concat(rows, ignore_index=True)
    # Coerce numeric intensity; drop NA/zero
    long_df['intensity'] = pd.to_numeric(long_df['intensity'], errors='coerce')
    long_df = long_df.dropna(subset=['intensity'])
    long_df = long_df[long_df['intensity'] > 0]
    Path(Path(args.out_csv).parent).mkdir(parents=True, exist_ok=True)
    long_df.to_csv(args.out_csv, index=False)
    print(f'Wrote {args.out_csv} with {len(long_df):,} rows')

if __name__ == '__main__':
    main()
