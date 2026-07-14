#!/usr/bin/env python3
"""
Rank NSP12 human interactors by integrating AP-MS evidence with AF3 complex interface metrics.

Inputs (defaults align with repo outputs):
- Interactors + structure table with iface metrics merged:
  projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ifaces.csv

Outputs:
- Ranked CSV: projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked.csv

Scoring (defaults, tunable via CLI):
  score = w_struct*S_struct + w_iptm*struct_iptm + w_ifaces*log1p(contact_pairs)
  where log base is normalized by log(max_norm) to keep term in ~[0,1].

Usage:
  python projects/unbiased_nsp12_20251015/scripts/24_rank_nsp12_interactors_with_structure.py \
    --in_csv projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ifaces.csv \
    --out_csv projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked.csv \
    --w_struct 0.6 --w_iptm 0.3 --w_ifaces 0.1 --ifaces_norm 200
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def to_float(v, default=None):
    try:
        return float(v)
    except Exception:
        return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in_csv', default='projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ifaces.csv')
    ap.add_argument('--out_csv', default='projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors_struct_ranked.csv')
    ap.add_argument('--w_struct', type=float, default=0.5)
    ap.add_argument('--w_iptm', type=float, default=0.25)
    ap.add_argument('--w_ifaces', type=float, default=0.1)
    ap.add_argument('--w_area', type=float, default=0.1)
    ap.add_argument('--w_bonds', type=float, default=0.05)
    ap.add_argument('--w_hydroph', type=float, default=0.0)
    ap.add_argument('--ifaces_norm', type=float, default=200.0)
    ap.add_argument('--area_norm', type=float, default=4000.0)
    args = ap.parse_args()

    rows = []
    with Path(args.in_csv).open('r', encoding='utf-8') as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        for row in reader:
            cp = to_float(row.get('iface_contact_pairs'), 0.0)
            iptm = to_float(row.get('struct_iptm'), 0.0)
            sstruct = to_float(row.get('S_struct'), 0.0)
            area = to_float(row.get('iface_area'), 0.0)
            hb = to_float(row.get('iface_hbonds'), 0.0)
            sb = to_float(row.get('iface_salt_bridges'), 0.0)
            hp = to_float(row.get('iface_hydroph_contacts'), 0.0)
            term_ifaces = 0.0
            if cp and args.ifaces_norm > 1:
                term_ifaces = math.log1p(cp) / math.log(args.ifaces_norm)
            term_area = 0.0
            if area and args.area_norm > 0:
                term_area = max(0.0, min(1.0, area / args.area_norm))
            term_bonds = 0.0
            if hb is not None or sb is not None:
                hb_n = (hb or 0.0) / 50.0
                sb_n = (sb or 0.0) / 20.0
                term_bonds = max(0.0, min(1.0, 0.6 * hb_n + 0.4 * sb_n))
            term_hyd = 0.0
            if hp is not None:
                term_hyd = max(0.0, min(1.0, (hp or 0.0) / 50.0))
            score = (
                args.w_struct * (sstruct or 0.0) +
                args.w_iptm * (iptm or 0.0) +
                args.w_ifaces * term_ifaces +
                args.w_area * term_area +
                args.w_bonds * term_bonds +
                args.w_hydroph * term_hyd
            )
            row['score_struct'] = f"{score:.6f}"
            rows.append(row)
    # sort desc
    rows.sort(key=lambda r: to_float(r.get('score_struct'), 0.0), reverse=True)
    # add rank
    for i, r in enumerate(rows, start=1):
        r['priority_rank'] = i
    out_cols = (reader.fieldnames or []) + ['score_struct', 'priority_rank']
    with Path(args.out_csv).open('w', encoding='utf-8', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=out_cols)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


if __name__ == '__main__':
    main()
