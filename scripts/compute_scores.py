#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apms_long', required=True)
    ap.add_argument('--control_stats', required=True)
    ap.add_argument('--out_csv', required=True)
    ap.add_argument('--z_thr', type=float, default=3.0)
    ap.add_argument('--log2fc_thr', type=float, default=1.0)
    ap.add_argument('--pseudo_percentile', type=int, default=5)
    args = ap.parse_args()

    df = pd.read_csv(args.apms_long)
    ctrl = pd.read_csv(args.control_stats)
    pseudo = np.percentile(df['intensity'], args.pseudo_percentile)
    df['log2_intensity'] = np.log2(df['intensity'] + pseudo)

    # Merge control stats by prey
    ctrl = ctrl.rename(columns={'prey_id': 'Protein IDs'})
    merged = df.merge(ctrl, on='Protein IDs', how='left')

    # Compute robust z and log2FC vs controls
    mad = merged['ctrl_MAD'].fillna(1.0)
    med = merged['ctrl_log2_median'].fillna(merged['log2_intensity'].median())
    merged['z_ctrl'] = (merged['log2_intensity'] - med) / (1.4826 * mad.replace(0, 1e-6))
    merged['log2FC_vs_ctrl'] = merged['log2_intensity'] - med

    # Compute bait-specific background across other baits in same strain (CompPASS-like)
    # For each (strain, prey), compute mean and std across baits; Z for the bait is deviation from others
    def comp_z(sub):
        # sub: rows for a given (strain, prey)
        vals = sub['log2_intensity'].values
        mean = vals.mean()
        std = vals.std(ddof=0) if len(vals) > 1 else 1.0
        return (sub['log2_intensity'] - mean) / (std if std>0 else 1.0)

    merged['strain_prey'] = merged['strain'].astype(str) + '|' + merged['Protein IDs'].astype(str)
    comp_z_series = merged.groupby('strain_prey', group_keys=False).apply(comp_z)
    merged['compZ_within_strain'] = comp_z_series.values

    # Robust deviation vs other baits (median-based), stable when only one bait has detection
    grp = merged.groupby(['strain','Protein IDs'])['log2_intensity']
    med_other = grp.transform('median')
    mad_other = grp.transform(lambda s: np.median(np.abs(s - np.median(s))))
    merged['delta_vs_other_baits'] = merged['log2_intensity'] - med_other
    merged['robust_compZ'] = merged['delta_vs_other_baits'] / (1.4826 * mad_other.replace(0, np.nan))
    merged['robust_compZ'] = merged['robust_compZ'].fillna(0.0)

    # Specificity: 1 - detection frequency across other baits in same strain
    det = merged.groupby(['strain', 'Protein IDs'])['log2_intensity'].apply(lambda s: s.notna().sum()).rename('n_detect_baits')
    merged = merged.merge(det.reset_index(), on=['strain','Protein IDs'], how='left')
    n_baits_per_strain = merged.groupby('strain')['bait'].nunique().to_dict()
    merged['n_baits'] = merged['strain'].map(n_baits_per_strain)
    merged['detect_freq'] = merged['n_detect_baits'] / merged['n_baits'].replace(0,1)
    merged['specificity'] = 1.0 - merged['detect_freq'].clip(0,1)

    # WD-like score
    merged['WD_score'] = merged['compZ_within_strain'] * merged['specificity']

    # MiST-like score: A * S * R; define R via cross-variant consistency of presence for this (bait, prey)
    # Compute presence matrix across strains
    presence = merged.assign(present=1).pivot_table(index=['bait','Protein IDs'], columns='strain', values='present', aggfunc='max').fillna(0)
    # Reproducibility: fraction of variants with presence
    rep = presence.mean(axis=1).rename('R')
    rep_df = rep.reset_index()
    merged = merged.merge(rep_df, on=['bait','Protein IDs'], how='left')
    # Abundance: scaled log2_intensity per bait
    merged['Abundance'] = merged.groupby(['strain','bait'])['log2_intensity'].transform(lambda s: (s - s.mean())/(s.std(ddof=0) if s.std(ddof=0)>0 else 1.0))
    merged['MiST_like'] = (merged['Abundance'].clip(-3,3) * merged['specificity'].clip(0,1) * merged['R'].clip(0,1)).clip(lower=0)

    # Promiscuity score: penalize frequent detection across baits and high control ratio
    merged['promiscuity_score'] = merged['detect_freq'].clip(0,1) * 0.6 + merged['ctrl_detect_ratio'].fillna(0).clip(0,1) * 0.4

    # Revised high/mid confidence flags (v2):
    # - High: strong enrichment vs control, not promiscuous, acceptable bait-specific deviation (median-based)
    merged['is_high_conf'] = (
        (merged['z_ctrl'] >= args.z_thr) &
        (merged['log2FC_vs_ctrl'] >= args.log2fc_thr) &
        (merged['promiscuity_score'] <= 0.3) &
        (merged['robust_compZ'] >= -0.5)
    )
    # - Mid: at least two moderate signals, limited promiscuity
    mid_signals = (
        (merged['z_ctrl'] >= 2.0).astype(int) +
        (merged['log2FC_vs_ctrl'] >= 0.75).astype(int) +
        (merged['robust_compZ'] >= 0.0).astype(int) +
        (merged['MiST_like'] >= 0.7).astype(int)
    )
    merged['is_mid_conf_v2'] = (mid_signals >= 2) & (merged['promiscuity_score'] <= 0.5)

    out = merged.drop(columns=['strain_prey'])
    Path(Path(args.out_csv).parent).mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f'Wrote {args.out_csv} with {len(out):,} rows')

if __name__ == '__main__':
    main()
