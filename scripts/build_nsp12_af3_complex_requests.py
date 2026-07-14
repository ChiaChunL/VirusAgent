#!/usr/bin/env python3
"""
Build AF3 complex JSON requests for NSP12 (viral) with top human interactors (AP-MS).

Inputs:
- NSP12 variant sequence (from data/按基因分类/按基因分类/NSP12.txt)
- Interactors CSV (e.g., outputs/tables/nsp12_interactors.csv)
- Host protein FASTA directories (projects/.../structure/sequences, sequences/)

Outputs:
- JSON requests under projects/unbiased_nsp12_20251015/structure/jobs/af3/requests

Usage example:
  python projects/unbiased_nsp12_20251015/scripts/23_build_nsp12_af3_complex_requests.py \
    --variant Prototype \
    --top_k 50
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_variant_sequences(txt_path: Path) -> Dict[str, str]:
    """Parse NSP12 variant sequences from a FASTA-like TXT.

    Expected header pattern: Variant-... on header lines; sequence lines starting with numbers.
    Returns mapping variant_name -> sequence (uppercase).
    """
    header_re = re.compile(r'^([A-Za-z][A-Za-z0-9_-]*)-')
    seq_line_re = re.compile(r'^\s*\d+')
    seqs: Dict[str, str] = {}
    cur: Optional[str] = None
    buf: List[str] = []
    with txt_path.open('r', encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            m = header_re.match(line)
            if m:
                if cur and buf:
                    seqs[cur] = ''.join(buf).upper()
                cur = m.group(1)
                buf = []
                continue
            if seq_line_re.match(line):
                part = re.sub(r'^\s*\d+\s*', '', line)
                part = re.sub(r'\s+', '', part)
                buf.append(part)
    if cur and buf:
        seqs[cur] = ''.join(buf).upper()
    # handle "与Prototype相同" markers by copying Prototype if present
    proto = None
    for k, v in seqs.items():
        if k.lower().startswith('prototype'):
            proto = v
            break
    if proto:
        for k in list(seqs.keys()):
            if '与Prototype相同' in k:
                seqs[k] = proto
    return seqs


def read_fasta_sequence(fp: Path) -> Optional[str]:
    try:
        txt = fp.read_text(encoding='utf-8', errors='ignore')
    except Exception:
        return None
    lines = [ln.strip() for ln in txt.splitlines() if ln.strip()]
    lines = [ln for ln in lines if not ln.startswith(';')]
    if not lines or not lines[0].startswith('>'):
        return None
    seq = ''.join(ln for ln in lines[1:] if not ln.startswith('>'))
    return re.sub(r'\s+', '', seq).upper()


def find_host_sequence(uniprot: str, search_dirs: List[Path]) -> Optional[str]:
    candidates = [d / f"{uniprot}.fasta" for d in search_dirs]
    for fp in candidates:
        if fp.exists():
            seq = read_fasta_sequence(fp)
            if seq:
                return seq
    return None


def load_interactors(csv_path: Path, top_k: int) -> List[str]:
    """Return list of UniProt IDs for top interactors.
    Accept column name 'prey_uniprot' or 'Protein IDs'.
    """
    ids: List[str] = []
    with csv_path.open('r', encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            uniprot = row.get('prey_uniprot') or row.get('Protein IDs')
            if not uniprot:
                continue
            ids.append(uniprot.strip())
            if len(ids) >= top_k:
                break
    return ids


def build_af3_complex_json(name: str, chain_a_seq: str, chain_b_seq: str) -> Dict:
    return {
        "name": name,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": chain_a_seq}},
            {"protein": {"id": ["B"], "sequence": chain_b_seq}},
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--viral_gene', default='NSP12')
    ap.add_argument('--variant', default='Prototype', help='Single variant (deprecated when --variants is used)')
    ap.add_argument('--variants', nargs='*', default=None, help='Multiple variants, e.g., Prototype Alpha Delta Omicron')
    ap.add_argument('--viral_sequences_txt', default='data/按基因分类/按基因分类/NSP12.txt')
    ap.add_argument('--interactors_csv', default='projects/unbiased_nsp12_20251015/outputs/tables/nsp12_interactors.csv')
    ap.add_argument('--host_seq_dirs', nargs='*', default=[
        'projects/unbiased_nsp12_20251015/structure/sequences',
        'sequences'
    ])
    ap.add_argument('--out_dir', default='projects/unbiased_nsp12_20251015/structure/jobs/af3/requests')
    ap.add_argument('--top_k', type=int, default=50)
    args = ap.parse_args()

    viral_txt = Path(args.viral_sequences_txt)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    host_dirs = [Path(p) for p in args.host_seq_dirs]

    # determine variants list
    variants = args.variants if args.variants else [args.variant]

    # load NSP12 variant sequences (map variant_prefix -> seq)
    seqs = parse_variant_sequences(viral_txt)
    # normalize common alias/typos, e.g., 'Protoype' -> 'Prototype'
    for k in list(seqs.keys()):
        kl = k.lower()
        if kl.startswith('proto') and 'prototype' not in seqs:
            seqs['Prototype'] = seqs[k]

    # load interactors list once
    uniprots = load_interactors(Path(args.interactors_csv), args.top_k)
    logger.info(f"Top interactors to build: {len(uniprots)} (top_k={args.top_k}) for variants={variants}")

    total_built = 0
    for var in variants:
        viral_seq = None
        for k, v in seqs.items():
            if k.startswith(var):
                viral_seq = v
                break
        if not viral_seq:
            logger.warning(f"Skip variant {var}: sequence not found in {viral_txt}")
            continue
        built, missing = 0, []
        for uid in uniprots:
            host_seq = find_host_sequence(uid, host_dirs)
            if not host_seq:
                missing.append(uid)
                continue
            name = f"{var}-{args.viral_gene}__{uid}"
            obj = build_af3_complex_json(name, viral_seq, host_seq)
            fp = out_dir / f"{name}.json"
            fp.write_text(json.dumps(obj, ensure_ascii=False, indent=2))
            built += 1
        total_built += built
        logger.info(f"Variant {var}: built {built} requests at {out_dir}; missing sequences: {len(missing)}")
        if missing:
            logger.warning(f"{var} missing host sequences (first 10): {missing[:10]}{'...' if len(missing)>10 else ''}")
    logger.info(f"All variants done. Total built: {total_built}")


if __name__ == '__main__':
    main()
