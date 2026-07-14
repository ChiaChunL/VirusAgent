#!/usr/bin/env python3
"""
生成“讲故事”式的可视化面板（Storyline Figures），用于从“为何选择NSP12”→
“无偏基线排名”→“Epsilon 不确定性”→“最小校准”的完整叙事。

输出图形（保存到 outputs/figs/）：
  - pca_tsne_highlight.png          （PCA/TSNE 高亮 NSP12）
  - apms_volcano.png                （AP‑MS 火山图，RIG‑I/PSMA1 高亮）
  - topk_unbiased_bar.png           （Top‑K 无偏排名条形图，RIG‑I 高亮）
  - epsilon_distribution.png        （ε 分布直方+KDE，含推荐线）
  - rank_dumbbell.png               （哑铃图：无偏 vs 校准 排名变化）
  - radar_RIGI_vs_Top1_vs_Neg.png   （雷达图：RIG‑I vs 无偏Top1 vs 低分对照）
  - nsp12_hotspot_lollipop.png      （NSP12 界面热点棒棒糖）
  - nsp12_hotspot_conservation.png  （NSP12 热点+保守性着色线图）

所依赖的表：
  - PCA:  af3_analysis/pca_analysis_results.csv
  - APMS: projects/unbiased_nsp12_20251015/outputs/tables/prey_scores.csv
  - 排名: projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked.csv
         projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked_calibrated.csv
  - ε:    projects/unbiased_nsp12_20251015/outputs/tables/epsilon_estimate.json（可选）
  - 热点: af3_analysis/data_summary/nsp12_interface_hotspots.csv
          af3_analysis/data_summary/nsp12_hotspots_conservation.csv

注意：本脚本仅生成图形，不改变任何排序或参数；与无偏基线完全解耦。
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[3]
FIG_DIR = ROOT / 'projects' / 'unbiased_nsp12_20251015' / 'outputs' / 'figs'
PCA_CSV = ROOT / 'af3_analysis' / 'pca_analysis_results.csv'
APMS_CSV = ROOT / 'projects' / 'unbiased_nsp12_20251015' / 'outputs' / 'tables' / 'prey_scores.csv'
RANK_CSV = ROOT / 'projects' / 'unbiased_nsp12_20251015' / 'outputs' / 'tables' / 'nsp12_interactors_struct_ranked.csv'
RANK_CAL_CSV = ROOT / 'projects' / 'unbiased_nsp12_20251015' / 'outputs' / 'tables' / 'nsp12_interactors_struct_ranked_calibrated.csv'
EPS_JSON = ROOT / 'projects' / 'unbiased_nsp12_20251015' / 'outputs' / 'tables' / 'epsilon_estimate.json'
HOTSPOT_CSV = ROOT / 'af3_analysis' / 'data_summary' / 'nsp12_interface_hotspots.csv'
CONSCSV = ROOT / 'af3_analysis' / 'data_summary' / 'nsp12_hotspots_conservation.csv'


def _to_float(v, d=0.0):
    try:
        return float(v)
    except Exception:
        return d


def pca_tsne_highlight(pca_csv: Path, out_png: Path):
    if not pca_csv.exists():
        return
    df = pd.read_csv(pca_csv)
    # 约定列：variant, gene, pca_component1, pca_component2
    df['is_nsp12'] = (df['gene'].astype(str).str.upper() == 'NSP12')
    sns.set_theme(context='paper', style='whitegrid', font_scale=1.1)
    plt.figure(figsize=(6,5))
    # 其他点灰色半透明，NSP12 高亮红色
    others = df[~df['is_nsp12']]
    nsp12 = df[df['is_nsp12']]
    plt.scatter(others['pca_component1'], others['pca_component2'], c='lightgray', alpha=0.5, s=20, label='Others')
    plt.scatter(nsp12['pca_component1'], nsp12['pca_component2'], c='#d62728', alpha=0.9, s=40, label='NSP12')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA scatter highlighting NSP12')
    plt.legend(frameon=False, loc='best')
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=600)
    plt.close()


def apms_volcano(apms_csv: Path, out_png: Path, highlight_ids: Optional[List[str]] = None):
    """基于 prey_scores.csv 生成更稳健的火山图：
    - 仅使用 NSP12 bait 且非对照行（若可得），否则使用全体非对照行。
    - 按 Protein IDs 聚合：x=median(log2FC_vs_ctrl)，y=-log10(1-specificity_mean)。
    - 显著性阈值：log2FC>1 且 specificity_mean>0.7（可通过脚本参数再调）。
    """
    if not apms_csv.exists():
        return
    usecols = ['Protein IDs','Protein names','log2FC_vs_ctrl','specificity','bait','is_control','detect_freq','MiST_like']
    df = pd.read_csv(apms_csv, usecols=[c for c in usecols if c in pd.read_csv(apms_csv, nrows=0).columns])
    # 仅 NSP12 & 非对照，如无 bait 列则只筛非对照
    if 'bait' in df.columns:
        sub = df[(df['bait'].astype(str).str.upper()=='NSP12') & (df.get('is_control', False)==False)].copy()
        if sub.empty:
            sub = df[df.get('is_control', False)==False].copy()
    else:
        sub = df.copy()
    # 聚合（去重）
    agg = sub.groupby('Protein IDs').agg(
        log2FC_med=('log2FC_vs_ctrl','median'),
        spec_mean=('specificity','mean'),
        mist_mean=('MiST_like','mean') if 'MiST_like' in sub.columns else ('specificity','mean'),
        n=('Protein IDs','count')
    ).reset_index()
    eps = 1e-6
    agg['x'] = agg['log2FC_med']
    agg['y'] = -np.log10(1.0 - agg['spec_mean'].clip(0, 0.999999) + eps)
    sig = (agg['x'] > 1.0) & (agg['spec_mean'] > 0.7)

    sns.set_theme(context='paper', style='whitegrid', font_scale=1.1)
    plt.figure(figsize=(7.2,5.2))
    plt.scatter(agg.loc[~sig, 'x'], agg.loc[~sig, 'y'], c='lightgray', alpha=0.6, s=14, label='Not significant')
    plt.scatter(agg.loc[sig, 'x'], agg.loc[sig, 'y'], c='#1f77b4', alpha=0.8, s=20, label='Significant')

    # 高亮指定蛋白（去重后仅标一次）
    if highlight_ids:
        for uid in highlight_ids:
            m = agg[agg['Protein IDs'] == uid]
            if not m.empty:
                x = float(m['x'].iloc[0]); y = float(m['y'].iloc[0])
                plt.scatter([x], [y], c='#d62728', s=40, label=uid)
                plt.text(x, y, uid, fontsize=9, color='#d62728')

    plt.axvline(1.0, color='gray', linestyle='--', linewidth=1)
    plt.axhline(-math.log10(1-0.7+eps), color='gray', linestyle='--', linewidth=1)
    plt.xlabel('Median log2 Fold Change vs control (per prey)')
    plt.ylabel('-log10(1 - mean specificity)')
    plt.title('AP-MS Volcano (NSP12 bait; per-prey aggregated)')
    plt.legend(frameon=False, loc='upper left')
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=600)
    plt.close()


def read_rank_tables() -> Tuple[pd.DataFrame, pd.DataFrame]:
    df_unb = pd.read_csv(RANK_CSV)
    df_cal = pd.read_csv(RANK_CAL_CSV) if RANK_CAL_CSV.exists() else None
    return df_unb, df_cal


def topk_unbiased_bar(df_unb: pd.DataFrame, out_png: Path, k: int = 30):
    d = df_unb.copy()
    d['score_struct'] = d['score_struct'].astype(float)
    d = d.sort_values('score_struct', ascending=False).head(k)
    # 标签优先 Protein IDs
    labels = d['Protein IDs'].fillna(d.get('prey_uniprot', '')).astype(str).tolist()
    scores = d['score_struct'].tolist()
    colors = ['#ff7f0e' if lab == 'O95786' else '#1f77b4' for lab in labels]
    plt.figure(figsize=(8,6))
    y = np.arange(len(labels))
    plt.barh(y, scores, color=colors, alpha=0.85)
    plt.yticks(y, labels)
    plt.gca().invert_yaxis()
    plt.xlabel('S_unbiased (score_struct)')
    plt.title('Top-K Unbiased Ranking (RIG-I highlighted)')
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=600)
    plt.close()


def epsilon_distribution(df_unb: pd.DataFrame, out_png: Path, n: int = 500, concentration: float = 300.0, q: float = 0.75):
    # 与 32_visualize_epsilon_distribution.py 一致的项
    def compute_score(row: pd.Series, w: np.ndarray, ifaces_norm=200.0, area_norm=4000.0) -> float:
        sstruct = _to_float(row.get('S_struct'), 0.0)
        iptm = _to_float(row.get('struct_iptm'), 0.0)
        cp = _to_float(row.get('iface_contact_pairs'), 0.0)
        area = _to_float(row.get('iface_area'), 0.0)
        hb = _to_float(row.get('iface_hbonds'), 0.0)
        sb = _to_float(row.get('iface_salt_bridges'), 0.0)
        hp = _to_float(row.get('iface_hydroph_contacts'), 0.0)
        term_ifaces = (np.log1p(cp) / np.log(200.0)) if (cp and 200.0 > 1) else 0.0
        term_area = max(0.0, min(1.0, area/4000.0)) if area else 0.0
        hb_n = (hb or 0.0)/50.0; sb_n = (sb or 0.0)/20.0
        term_bonds = max(0.0, min(1.0, 0.6*hb_n + 0.4*sb_n))
        term_hyd = max(0.0, min(1.0, (hp or 0.0)/50.0))
        terms = np.array([sstruct, iptm, term_ifaces, term_area, term_bonds, term_hyd], dtype=float)
        return float(np.dot(w, terms))

    base = np.array([0.5, 0.25, 0.1, 0.1, 0.05, 0.0], dtype=float)
    base = base/base.sum()
    alpha = np.maximum(base * concentration, 1e-6)
    scores = df_unb.copy()
    gaps: List[float] = []
    for _ in range(n):
        w = np.random.dirichlet(alpha)
        vals = [compute_score(r, w) for _, r in scores.iterrows()]
        vals.sort(reverse=True)
        if len(vals) >= 2:
            gaps.append(vals[0]-vals[1])
    if not gaps:
        return
    gaps = np.array(gaps)
    eps_est = float(np.quantile(gaps, q))
    eps_rec = None
    if EPS_JSON.exists():
        try:
            eps_rec = json.loads(EPS_JSON.read_text()).get('recommended_epsilon')
        except Exception:
            eps_rec = None
    sns.set_theme(context='paper', style='whitegrid', font_scale=1.1)
    plt.figure(figsize=(6,4))
    sns.histplot(gaps, bins=30, kde=True, color='#1f77b4', alpha=0.7)
    plt.axvline(eps_est, color='#2ca02c', linestyle='--', label=f'Estimated q={q:.2f}: {eps_est:.4f}')
    if eps_rec is not None:
        plt.axvline(eps_rec, color='#d62728', linestyle='-.', label=f'Recommended: {eps_rec:.4f}')
    plt.xlabel('Top-2 score gap under weight perturbations')
    plt.ylabel('Frequency')
    plt.title('Epsilon (equivalence margin) distribution')
    plt.legend(frameon=False)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=600)
    plt.close()


def rank_dumbbell(df_unb: pd.DataFrame, df_cal: Optional[pd.DataFrame], out_png: Path, top_n: int = 15):
    if df_cal is None or df_cal.empty:
        return
    pre_rank = {str(r.get('Protein IDs') or r.get('prey_uniprot')): int(r.get('priority_rank') or 0) for _, r in df_unb.iterrows()}
    post_rank = {str(r.get('Protein IDs') or r.get('prey_uniprot')): int(r.get('priority_rank_cal') or 0) for _, r in df_cal.iterrows()}
    ids = [i for i in pre_rank.keys() if i in post_rank]
    # 仅显示排名变化者 + 无偏Top
    ids = sorted(ids, key=lambda x: pre_rank[x])[:top_n]
    pre_vals = np.array([pre_rank[i] for i in ids])
    post_vals = np.array([post_rank[i] for i in ids])
    y = np.arange(len(ids))
    sns.set_theme(context='paper', style='whitegrid', font_scale=1.1)
    plt.figure(figsize=(8,5))
    for i in range(len(ids)):
        c = '#d62728' if ids[i] == 'O95786' else '#7f7f7f'
        plt.plot([pre_vals[i], post_vals[i]], [y[i], y[i]], '-', color=c, alpha=0.7)
        plt.scatter(pre_vals[i], y[i], c='#1f77b4', s=30, label='Unbiased' if i==0 else None)
        plt.scatter(post_vals[i], y[i], c='#2ca02c', s=30, label='Calibrated' if i==0 else None)
    plt.yticks(y, ids)
    plt.gca().invert_xaxis()
    plt.xlabel('Rank (smaller is better)')
    plt.title('Rank change after minimal prior calibration')
    plt.legend(frameon=False, loc='upper right')
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=600)
    plt.close()


def _norm_terms_row(r: pd.Series) -> Dict[str, float]:
    sstruct = _to_float(r.get('S_struct'), 0.0)
    iptm = _to_float(r.get('struct_iptm'), 0.0)
    cp = _to_float(r.get('iface_contact_pairs'), 0.0)
    area = _to_float(r.get('iface_area'), 0.0)
    hb = _to_float(r.get('iface_hbonds'), 0.0)
    sb = _to_float(r.get('iface_salt_bridges'), 0.0)
    plddt = _to_float(r.get('iface_mean_plddt'), 0.0)
    contacts_norm = (math.log1p(cp) / math.log(200.0)) if (cp and 200.0 > 1) else 0.0
    area_norm = max(0.0, min(1.0, area/4000.0)) if area else 0.0
    bonds_norm = max(0.0, min(1.0, 0.6*(hb/50.0) + 0.4*(sb/20.0)))
    plddt_norm = max(0.0, min(1.0, plddt/100.0))
    return {
        'S_struct': sstruct,
        'ipTM': iptm,
        'contacts_norm': contacts_norm,
        'area_norm': area_norm,
        'bonds_norm': bonds_norm,
        'iface_pLDDT': plddt_norm,
    }


def radar_plot(df_unb: pd.DataFrame, out_png: Path):
    d = df_unb.copy()
    d['score_struct'] = d['score_struct'].astype(float)
    d = d.sort_values('score_struct', ascending=False)
    # 挑选：RIG‑I(O95786)，无偏Top1（若非RIG‑I），以及一个低分样本
    rig = d[d['Protein IDs'].astype(str) == 'O95786'].head(1)
    top1 = d.head(1)
    if not rig.empty and (top1['Protein IDs'].iloc[0] == 'O95786'):
        top1 = d.iloc[[1]]  # 若 RIG‑I 为Top1，改取第二名
    neg = d.tail(1)

    candidates = []
    labels = []
    for sub, name in [(rig, 'RIG-I (O95786)'), (top1, f"Top1 ({top1['Protein IDs'].iloc[0]})"), (neg, f"Low ({neg['Protein IDs'].iloc[0]})")]:
        if sub is not None and not sub.empty:
            labels.append(name)
            candidates.append(_norm_terms_row(sub.iloc[0]))

    if not candidates:
        return
    metrics = ['S_struct', 'ipTM', 'contacts_norm', 'area_norm', 'bonds_norm', 'iface_pLDDT']
    theta = np.linspace(0, 2*np.pi, len(metrics), endpoint=False)
    sns.set_theme(context='paper', style='whitegrid', font_scale=1.0)
    plt.figure(figsize=(6.2,6.2))
    ax = plt.subplot(111, polar=True)
    for i, c in enumerate(candidates):
        vals = [c[m] for m in metrics]
        vals += vals[:1]
        tt = np.concatenate([theta, theta[:1]])
        color = ['#d62728', '#1f77b4', '#7f7f7f'][i % 3]
        ax.plot(tt, vals, color=color, linewidth=2, label=labels[i])
        ax.fill(tt, vals, color=color, alpha=0.1)
    ax.set_thetagrids(theta * 180/np.pi, metrics)
    ax.set_ylim(0, 1)
    plt.title('Feature radar: RIG-I vs Unbiased Top1 vs Low')
    plt.legend(loc='upper right', bbox_to_anchor=(1.35, 1.1), frameon=False)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close()


def hotspot_plots(hot_csv: Path, cons_csv: Path, out1: Path, out2: Path, top_k: int = 25):
    if hot_csv.exists():
        hot = pd.read_csv(hot_csv)
        items = hot[['resA_idx', 'contact_count']].dropna()
        items = items.sort_values('contact_count', ascending=False).head(top_k)
        plt.figure(figsize=(9,3.6))
        # 兼容旧版 matplotlib：不使用 use_line_collection 参数
        plt.stem(items['resA_idx'], items['contact_count'], basefmt=' ')
        plt.xlabel('NSP12 residue index')
        plt.ylabel('Contact count')
        plt.title('NSP12 interface hotspots (Top-25)')
        plt.tight_layout()
        out1.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out1, dpi=600)
        plt.close()
    if cons_csv.exists():
        cons = pd.read_csv(cons_csv)
        cons = cons[['position', 'contact_count', 'conservation_score']].dropna()
        cons = cons.sort_values('position')
        x = cons['position'].values
        y = cons['contact_count'].values
        c = cons['conservation_score'].values
        plt.figure(figsize=(10,3.6))
        sc = plt.scatter(x, y, c=c, cmap='coolwarm', s=18)
        plt.colorbar(sc, label='Conservation score')
        plt.xlabel('NSP12 residue position')
        plt.ylabel('Contact count')
        plt.title('NSP12 interface hotspots with conservation coloring')
        plt.tight_layout()
        out2.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out2, dpi=600)
        plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fig_dir', default=str(FIG_DIR))
    ap.add_argument('--topk', type=int, default=30)
    ap.add_argument('--eps_draws', type=int, default=500)
    ap.add_argument('--eps_concentration', type=float, default=300.0)
    ap.add_argument('--eps_q', type=float, default=0.75)
    args = ap.parse_args()

    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 1) PCA/TSNE 高亮 NSP12
    pca_tsne_highlight(PCA_CSV, fig_dir / 'pca_tsne_highlight.png')

    # 2) AP‑MS 火山图：高亮 RIG‑I + 无偏Top1（若 Top1 就是 RIG‑I，则高亮 Top2）
    top_ids = []
    try:
        dtmp = df_unb.sort_values('score_struct', ascending=False).reset_index(drop=True)
        if not dtmp.empty:
            t1 = str(dtmp.loc[0, 'Protein IDs'])
            top_ids.append(t1)
            if t1 == 'O95786' and len(dtmp) > 1:
                top_ids.append(str(dtmp.loc[1, 'Protein IDs']))
        if 'O95786' not in top_ids:
            top_ids.append('O95786')
    except Exception:
        top_ids = ['O95786']
    apms_volcano(APMS_CSV, fig_dir / 'apms_volcano.png', highlight_ids=top_ids)

    # 3) Top‑K 无偏排名条形图
    df_unb, df_cal = read_rank_tables()
    topk_unbiased_bar(df_unb, fig_dir / 'topk_unbiased_bar.png', k=args.topk)

    # 4) ε 分布
    epsilon_distribution(df_unb, fig_dir / 'epsilon_distribution.png', n=args.eps_draws, concentration=args.eps_concentration, q=args.eps_q)

    # 5) 哑铃图（排名变化）
    rank_dumbbell(df_unb, df_cal, fig_dir / 'rank_dumbbell.png', top_n=15)

    # 6) 雷达图（RIG‑I vs 无偏Top1 vs 低分）
    radar_plot(df_unb, fig_dir / 'radar_RIGI_vs_Top1_vs_Neg.png')

    # 7) NSP12 热点
    hotspot_plots(HOTSPOT_CSV, CONSCSV, fig_dir / 'nsp12_hotspot_lollipop.png', fig_dir / 'nsp12_hotspot_conservation.png', top_k=25)

    print('Storyline figures written to', fig_dir)


if __name__ == '__main__':
    main()
