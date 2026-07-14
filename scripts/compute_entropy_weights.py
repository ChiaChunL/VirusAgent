#!/usr/bin/env python3
"""
Compute objective feature weights using entropy weight method (EWM).
This is feature-level (criteria) weighting for multi-criteria ranking (e.g., TOPSIS).
Unsupervised, no labels, unbiased by protein identity.
"""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd


def entropy_weights(df: pd.DataFrame, feature_cols):
    X = df[feature_cols].values.astype(float)
    # Min-max normalize to [0,1]; if constant, set to small epsilon
    minv = X.min(axis=0)
    maxv = X.max(axis=0)
    rng = np.where((maxv - minv) == 0, 1.0, (maxv - minv))
    Xn = (X - minv) / rng
    # Avoid zeros for log; normalize per-feature so rows sum to 1
    eps = 1e-12
    P = Xn / (Xn.sum(axis=0) + eps)
    P = np.where(P <= 0, eps, P)
    m = Xn.shape[0]
    k = 1.0 / np.log(m + 1e-12)
    E = -k * (P * np.log(P)).sum(axis=0)  # entropy per feature, in [0,1]
    d = 1.0 - E  # diversity
    if d.sum() == 0:
        w = np.ones_like(d) / len(d)
    else:
        w = d / d.sum()
    return {col: float(w[i]) for i, col in enumerate(feature_cols)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features_csv', required=True)
    ap.add_argument('--out_json', required=True)
    args = ap.parse_args()

    feats = pd.read_csv(args.features_csv)
    id_cols = ['strain','bait','protein']
    base_cols = ['n_high','mean_WD','max_z','mean_MiST','deg_centrality','eig_centrality','score_pca1','score_iforest']
    feature_cols = [c for c in base_cols if c in feats.columns]
    if not feature_cols:
        raise SystemExit('No feature columns available for weighting.')

    w = entropy_weights(feats, feature_cols)
    Path(Path(args.out_json).parent).mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps({'feature_weights': w, 'method': 'entropy', 'feature_set': feature_cols}, indent=2))
    print('Saved weights to', args.out_json)

if __name__ == '__main__':
    main()

