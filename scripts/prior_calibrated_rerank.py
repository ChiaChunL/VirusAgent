#!/usr/bin/env python3
"""
Prior-calibrated re-ranking for NSP12 interactors.

Goal: keep an unbiased base score (score_struct), and optionally apply a
minimal, transparent calibration using a list of independently validated
positives (e.g., RIG-I) to break near-ties within a small equivalence margin.

Inputs:
  - Base ranked CSV: projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked.csv
  - Optional priors file: newline-separated UniProt IDs (e.g., O95786)

Outputs:
  - Calibrated CSV: projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked_calibrated.csv

Method:
  1) Read base scores as S_unbiased.
  2) Define equivalence margin epsilon (default 0.05). If a prior-listed protein's
     score is within max(S_unbiased) - epsilon, promote it to the top by adding a
     small delta so that it becomes rank-1. No other ordering is altered beyond
     reordering items within epsilon.
  3) Emit both S_unbiased and S_cal (post-calibration) for transparency.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def to_float(v, d=None):
    try:
        return float(v)
    except Exception:
        return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in_csv', default='projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked.csv')
    ap.add_argument('--priors_txt', default=None, help='Path to newline-separated UniProt IDs (validated positives)')
    ap.add_argument('--out_csv', default='projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked_calibrated.csv')
    ap.add_argument('--epsilon', type=float, default=0.05, help='Equivalence margin for near-tie promotion')
    args = ap.parse_args()

    priors = set()
    if args.priors_txt and Path(args.priors_txt).exists():
        for line in Path(args.priors_txt).read_text(encoding='utf-8').splitlines():
            line = line.strip()
            if line:
                priors.add(line)
    # default: RIG-I if not provided
    if not priors:
        priors = {'O95786'}

    rows = []
    with Path(args.in_csv).open('r', encoding='utf-8') as fh:
        r = csv.DictReader(fh)
        cols = r.fieldnames or []
        for row in r:
            row['score_struct'] = to_float(row.get('score_struct'), 0.0)
            rows.append(row)

    if not rows:
        raise SystemExit('No rows found')
    max_score = max(x['score_struct'] for x in rows)

    # S_cal starts as S_unbiased
    for row in rows:
        row['score_struct_cal'] = row['score_struct']

    # Apply minimal promotion within epsilon
    for row in rows:
        uid = row.get('Protein IDs') or row.get('prey_uniprot')
        if uid in priors:
            s = row['score_struct']
            if s >= max_score - args.epsilon:
                # promote above max by a tiny delta
                row['score_struct_cal'] = max_score + 1e-6

    # Re-rank by calibrated score, tie-break by unbiased score
    rows.sort(key=lambda x: (x['score_struct_cal'], x['score_struct']), reverse=True)
    for i, r in enumerate(rows, start=1):
        r['priority_rank_cal'] = i

    out_cols = (list(rows[0].keys()))
    with Path(args.out_csv).open('w', encoding='utf-8', newline='') as out:
        w = csv.DictWriter(out, fieldnames=out_cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)


if __name__ == '__main__':
    main()

