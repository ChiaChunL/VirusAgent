#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import json

try:
    import networkx as nx
except Exception:
    nx = None


def load_structure_features(highdim_dir: Path) -> pd.DataFrame:
    meta_path = Path(highdim_dir, 'metadata.json')
    feats_path = Path(highdim_dir, 'raw_features.npy')
    if not meta_path.exists() or not feats_path.exists():
        # Graceful fallback: return empty DF with required id columns
        return pd.DataFrame(columns=['strain','protein'])
    meta = json.loads(meta_path.read_text())
    names = meta.get('structure_names', [])
    feats = np.load(feats_path)
    feat_names = meta.get('feature_names', [f'feat_{i}' for i in range(int(feats.shape[1]))])
    rows = []
    for n, f in zip(names, feats):
        if isinstance(n, str) and (n.startswith('seed-') or n.startswith('sample-')):
            continue
        label = str(n).replace('Protoype_', 'Prototype_')
        row = {'label': label}
        for i, val in enumerate(f):
            key = feat_names[i] if i < len(feat_names) else f'feat_{i}'
            row[key] = float(val)
        rows.append(row)
    if not rows:
        return pd.DataFrame(columns=['strain','protein'])
    df = pd.DataFrame(rows)
    def split_label(s: str):
        if '_' in s:
            strain, protein = s.split('_', 1)
        elif '-' in s:
            strain, protein = s.split('-', 1)
        else:
            strain, protein = 'Unknown', s
        return strain, protein
    lab = df['label'].apply(split_label)
    df['strain'] = lab.apply(lambda x: x[0])
    df['protein'] = lab.apply(lambda x: x[1])
    agg = df.drop(columns=['label']).groupby(['strain','protein']).mean(numeric_only=True).reset_index()
    return agg


def build_network_features(scores: pd.DataFrame) -> pd.DataFrame:
    if nx is None:
        # Minimal degenerate features if networkx not available
        return scores.groupby(['strain','bait']).agg(
            n_prey=('Protein IDs','nunique'),
            n_high=('is_high_conf','sum'),
            mean_z=('z_ctrl','mean'),
            max_z=('z_ctrl','max'),
            mean_WD=('WD_score','mean'),
        ).reset_index()
    df = scores.copy()
    feats = []
    for strain, sub in df.groupby('strain'):
        G = nx.Graph()
        # add bait nodes and prey nodes
        baits = sub['bait'].unique().tolist()
        preys = sub['Protein IDs'].unique().tolist()
        G.add_nodes_from([(f'B:{b}', {'type': 'bait'}) for b in baits])
        G.add_nodes_from([(f'P:{p}', {'type': 'prey'}) for p in preys])
        # edges weighted by MiST-like or WD
        for _, r in sub.iterrows():
            G.add_edge(f'B:{r.bait}', f'P:{r["Protein IDs"]}', weight=max(0.0, float(r['MiST_like'])))
        # compute centrality for bait nodes
        deg = nx.degree_centrality(G)
        eig = nx.eigenvector_centrality_numpy(G, weight='weight') if G.number_of_edges() > 0 else {n:0 for n in G.nodes}
        for b in baits:
            node = f'B:{b}'
            feats.append({
                'strain': strain,
                'bait': b,
                'deg_centrality': deg.get(node, 0.0),
                'eig_centrality': eig.get(node, 0.0),
            })
    fdf = pd.DataFrame(feats)
    agg = scores.groupby(['strain','bait']).agg(
        n_prey=('Protein IDs','nunique'),
        n_high=('is_high_conf','sum'),
        mean_z=('z_ctrl','mean'),
        max_z=('z_ctrl','max'),
        mean_WD=('WD_score','mean'),
        mean_MiST=('MiST_like','mean'),
        max_MiST=('MiST_like','max'),
    ).reset_index()
    out = agg.merge(fdf, on=['strain','bait'], how='left')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scores_csv', required=True)
    ap.add_argument('--highdim_dir', required=True)
    ap.add_argument('--out_csv', required=True)
    args = ap.parse_args()

    scores = pd.read_csv(args.scores_csv)
    # Aggregate to gene-level per (strain, bait)
    net_feats = build_network_features(scores)

    # Map bait names to protein labels for highdim merge
    net_feats['protein'] = net_feats['bait']
    struct = load_structure_features(Path(args.highdim_dir))
    merged = net_feats.merge(struct, on=['strain','protein'], how='left')
    Path(Path(args.out_csv).parent).mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out_csv, index=False)
    print(f'Wrote features matrix: {args.out_csv} with {len(merged):,} rows')

if __name__ == '__main__':
    main()
