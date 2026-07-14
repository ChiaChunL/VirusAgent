#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest


def rank_aggregate(df: pd.DataFrame, cols):
    ranks = []
    for c in cols:
        asc = False  # higher is better
        r = df[c].rank(ascending=not asc, method='average')
        ranks.append(r)
    R = sum(ranks) / len(cols)
    return R


def topsis(df: pd.DataFrame, cols, weights=None):
    X = df[cols].values.astype(float)
    # vector normalization
    norm = np.linalg.norm(X, axis=0)
    norm[norm == 0] = 1.0
    Xn = X / norm
    w = np.ones(len(cols)) if weights is None else np.array(weights)
    w = w / w.sum()
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
    ap.add_argument('--out_csv', required=True)
    ap.add_argument('--eval_txt', required=True)
    args = ap.parse_args()

    feats = pd.read_csv(args.features_csv)
    # Exclude controls (bait/protein labeled as CON)
    feats = feats[feats['protein'].str.upper() != 'CON'].copy()
    # Select numeric columns for PCA and TOPSIS
    id_cols = ['strain','bait','protein']
    num_cols = [c for c in feats.columns if c not in id_cols and feats[c].dtype != 'object']
    X = feats[num_cols].fillna(feats[num_cols].median())

    # Standardize
    Xs = (X - X.mean()) / X.std(ddof=0).replace(0,1)

    # PCA with orientation: make PC1 positively correlate with beneficial feature if available
    pca = PCA(n_components=min(3, Xs.shape[1]))
    PCs = pca.fit_transform(Xs.values)
    score_pca1 = PCs[:,0]
    orient_ref = None
    for ref in ['n_high','mean_WD','max_z','mean_MiST']:
        if ref in feats.columns:
            orient_ref = feats[ref].values
            break
    if orient_ref is not None:
        corr = np.corrcoef(score_pca1, orient_ref)[0,1]
        if np.isnan(corr):
            corr = 1.0
        score_pca1 = np.sign(corr if corr != 0 else 1.0) * score_pca1
    feats['score_pca1'] = score_pca1

    # Isolation Forest (unsupervised anomaly)
    iforest = IsolationForest(random_state=42, contamination='auto')
    iforest.fit(Xs.values)
    feats['score_iforest'] = -iforest.score_samples(Xs.values)  # higher=more anomalous

    # Aggregated ranking
    base_cols = ['n_high','mean_WD','max_z','mean_MiST','deg_centrality','eig_centrality','score_pca1','score_iforest']
    base_cols = [c for c in base_cols if c in feats.columns]
    feats['score_rankagg'] = rank_aggregate(feats, base_cols)

    # TOPSIS with uniform weights
    feats['score_topsis'] = topsis(feats, base_cols)

    # Final composite (simple average of normalized scores)
    agg = feats[base_cols].apply(lambda s: (s - s.mean())/(s.std(ddof=0) if s.std(ddof=0)>0 else 1.0))
    feats['score_composite'] = agg.mean(axis=1)

    # Rank within-strain and global
    feats['rank_global'] = feats['score_composite'].rank(ascending=False, method='min')
    feats['rank_within_strain'] = feats.groupby('strain')['score_composite'].rank(ascending=False, method='min')

    # Save
    Path(Path(args.out_csv).parent).mkdir(parents=True, exist_ok=True)
    feats.to_csv(args.out_csv, index=False)

    # Evaluation: locate NSP12 positions without any special handling
    report = []
    sub = feats[feats['protein'].str.upper()=='NSP12']
    if len(sub) > 0:
        for _, r in sub.iterrows():
            report.append(f"NSP12 @ {r['strain']}: global_rank={int(r['rank_global'])}, within_strain_rank={int(r['rank_within_strain'])}, composite={r['score_composite']:.3f}")
    else:
        report.append('NSP12 not found in features table.')
    Path(Path(args.eval_txt).parent).mkdir(parents=True, exist_ok=True)
    Path(args.eval_txt).write_text('\n'.join(report))
    print('\n'.join(report))

if __name__ == '__main__':
    main()
