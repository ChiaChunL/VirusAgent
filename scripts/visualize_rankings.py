#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style='whitegrid', context='talk', font_scale=0.9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ranking_csv', required=True)
    ap.add_argument('--figs_dir', required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.ranking_csv)
    figs_dir = Path(args.figs_dir)
    figs_dir.mkdir(parents=True, exist_ok=True)

    # Global top 10 bar
    top = df.sort_values('score_composite', ascending=False).head(10)
    plt.figure(figsize=(10,4))
    plt.bar(top['protein'] + '_' + top['strain'], top['score_composite'])
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(figs_dir / 'global_top10.png', dpi=300); plt.savefig(figs_dir / 'global_top10.pdf')

    # NSP12 profile if present
    nsp = df[df['protein'].str.upper()=='NSP12']
    if len(nsp) > 0:
        plt.figure(figsize=(6,4))
        plt.bar(nsp['strain'], nsp['score_composite'])
        plt.ylabel('Composite score')
        plt.title('NSP12 across variants')
        plt.tight_layout()
        plt.savefig(figs_dir / 'nsp12_by_variant.png', dpi=300); plt.savefig(figs_dir / 'nsp12_by_variant.pdf')

if __name__ == '__main__':
    main()
