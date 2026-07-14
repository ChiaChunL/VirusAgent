#!/usr/bin/env python3
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def topsis(df: pd.DataFrame, cols, weights):
    X = df[cols].values.astype(float)
    norm = np.linalg.norm(X, axis=0)
    norm[norm == 0] = 1.0
    Xn = X / norm
    w = np.array([weights.get(c, 1.0) for c in cols], dtype=float)
    w = w / (w.sum() if w.sum() > 0 else 1.0)
    V = Xn * w
    ideal = V.max(axis=0)
    nadir = V.min(axis=0)
    d_pos = np.linalg.norm(V - ideal, axis=1)
    d_neg = np.linalg.norm(V - nadir, axis=1)
    score = d_neg / (d_pos + d_neg + 1e-9)
    return score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features_csv', required=True)
    ap.add_argument('--weights_json', required=True)
    ap.add_argument('--out_csv', required=True)
    ap.add_argument('--eval_txt', required=True)
    args = ap.parse_args()

    feats = pd.read_csv(args.features_csv)
    with open(args.weights_json, 'r') as f:
        wobj = json.load(f)
    weights = wobj['feature_weights']
    features_used = wobj['feature_set']

    # Exclude controls
    feats = feats[feats['protein'].str.upper() != 'CON'].copy()
    id_cols = ['strain','bait','protein']
    num_cols = [c for c in feats.columns if c not in id_cols and feats[c].dtype != 'object']
    X = feats[num_cols].fillna(feats[num_cols].median())
    Xs = (X - X.mean()) / X.std(ddof=0).replace(0,1)

    # PCA for context (not used directly in weighted score)
    pca = PCA(n_components=min(3, Xs.shape[1]))
    PCs = pca.fit_transform(Xs.values)
    feats['score_pca1_recomp'] = PCs[:,0]

    # Weighted TOPSIS
    cols = [c for c in features_used if c in feats.columns]
    feats['score_topsis_w'] = topsis(feats, cols, weights)

    # Weighted composite: normalized features weighted by w
    Z = feats[cols].apply(lambda s: (s - s.mean())/(s.std(ddof=0) if s.std(ddof=0)>0 else 1.0))
    wvec = np.array([weights.get(c,1.0) for c in cols])
    wvec = wvec / (wvec.sum() if wvec.sum()>0 else 1.0)
    feats['score_composite_w'] = (Z.values * wvec).sum(axis=1)

    # Ranks
    feats['rank_global_w'] = feats['score_composite_w'].rank(ascending=False, method='min')
    feats['rank_within_strain_w'] = feats.groupby('strain')['score_composite_w'].rank(ascending=False, method='min')

    Path(Path(args.out_csv).parent).mkdir(parents=True, exist_ok=True)
    feats.to_csv(args.out_csv, index=False)

    sub = feats[feats['protein'].str.upper()=='NSP12']
    report = []
    if len(sub) > 0:
        for _, r in sub.iterrows():
            report.append(f"[W] NSP12 @ {r['strain']}: global_rank={int(r['rank_global_w'])}, within_strain_rank={int(r['rank_within_strain_w'])}, composite_w={r['score_composite_w']:.3f}")
    else:
        report.append('[W] NSP12 not found in features table.')
    Path(Path(args.eval_txt).parent).mkdir(parents=True, exist_ok=True)
    Path(args.eval_txt).write_text('\n'.join(report))
    print('\n'.join(report))

if __name__ == '__main__':
    main()

