#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

sns.set_theme(style='whitegrid', context='talk', font_scale=0.9)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prey_scores_csv', required=True)
    ap.add_argument('--figs_dir', required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.prey_scores_csv)
    figs_dir = Path(args.figs_dir)

    # --- Volcano Plot ---
    df['neg_log10_z'] = df['z_ctrl']
    df_plot = df.drop_duplicates(subset=['Protein names']).copy()

    plt.figure(figsize=(12, 8))
    sns.scatterplot(data=df_plot, x='log2FC_vs_ctrl', y='neg_log10_z', alpha=0.5, edgecolor='none')

    # Highlight significant points
    significant = df_plot[(df_plot['log2FC_vs_ctrl'].abs() > 1) & (df_plot['neg_log10_z'] > 3)]
    sns.scatterplot(data=significant, x='log2FC_vs_ctrl', y='neg_log10_z', color='red', s=50, edgecolor='black', linewidth=0.5)

    # Annotate top proteins
    top_n = 10
    df_plot['rank_score'] = df_plot['log2FC_vs_ctrl'].abs() * df_plot['neg_log10_z']
    top_proteins = df_plot.nlargest(top_n, 'rank_score')

    # Add RIG-I if not in top N
    rig_i_protein = df_plot[df_plot['Protein names'].str.contains('DDX58', na=False)] # DDX58 is the gene name for RIG-I
    if not rig_i_protein.empty:
        top_proteins = pd.concat([top_proteins, rig_i_protein]).drop_duplicates()

    texts = []
    for i, row in top_proteins.iterrows():
        texts.append(plt.text(row['log2FC_vs_ctrl'], row['neg_log10_z'], row['Protein names'], fontsize=9))

    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))

    plt.title('Volcano Plot of NSP12 Interactors', fontsize=16)
    plt.xlabel('Log2 Fold Change vs Control', fontsize=12)
    plt.ylabel('Significance (Z-score vs Control)', fontsize=12)
    plt.axhline(3, ls='--', color='gray')
    plt.axvline(1, ls='--', color='gray')
    plt.axvline(-1, ls='--', color='gray')
    plt.tight_layout()
    plt.savefig(figs_dir / 'volcano_plot.png', dpi=300)
    plt.savefig(figs_dir / 'volcano_plot.pdf')

    print(f"Volcano plot saved to {figs_dir}")

if __name__ == '__main__':
    main()
