#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def robust_stats(x: pd.Series):
    x = x.dropna()
    if len(x) == 0:
        return np.nan, np.nan
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        mad = 1e-6
    return med, mad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apms_long', required=True)
    ap.add_argument('--out_csv', required=True)
    ap.add_argument('--pseudo_percentile', type=int, default=5)
    args = ap.parse_args()

    df = pd.read_csv(args.apms_long)
    # Pseudo for log2
    pseudo = np.percentile(df['intensity'], args.pseudo_percentile)
    df['log2_intensity'] = np.log2(df['intensity'] + pseudo)

    # Build control stats per prey across controls (aggregate across all control files)
    ctrl = df[df['is_control']]
    detect = (ctrl.groupby('Protein IDs')['log2_intensity']
                 .agg(['count', 'median'])
                 .rename(columns={'count': 'ctrl_detect_count', 'median': 'ctrl_log2_median'}))
    # MAD per prey across controls
    mad_series = ctrl.groupby('Protein IDs')['log2_intensity'].apply(lambda s: np.median(np.abs(s - np.median(s))) if len(s)>0 else np.nan)
    detect['ctrl_MAD'] = mad_series

    # Detection ratio: count / number of control files
    n_controls = ctrl[['strain','bait','source_file']].drop_duplicates().shape[0]
    detect['ctrl_detect_ratio'] = detect['ctrl_detect_count'] / max(1, n_controls)

    out = detect.reset_index().rename(columns={'Protein IDs': 'prey_id'})
    Path(Path(args.out_csv).parent).mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f'Wrote {args.out_csv} with {len(out):,} preys; controls={n_controls}')

if __name__ == '__main__':
    main()
