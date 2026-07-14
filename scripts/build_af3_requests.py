#!/usr/bin/env python3
"""
Build AlphaFold3 JSON requests from gene-variant text inputs.

Source data: data/按基因分类/按基因分类/*.txt
Output JSON: af3/input/<Variant_Gene>_input.json

Usage:
  python final/scripts/01_build_af3_requests_from_genes.py \
    --data_dir data/按基因分类/按基因分类 \
    --out_dir af3/input

Notes:
  - This script only prepares AF3 JSON requests; submission is separate.
  - Parsing mirrors monitor_af3.py to handle "与Prototype相同" cases.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path
from typing import Dict

# Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("build_af3_requests")

HEADER_PATTERN = re.compile(r'^([A-Za-z][A-Za-z0-9_-]*)-([A-Za-z0-9_-]+)')
SEQUENCE_LINE_PATTERN = re.compile(r'^\s*\d+')


def parse_gene_file(path: Path) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    same_as_proto = []
    current_variant = None
    parts = []

    with path.open('r', encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            m = HEADER_PATTERN.match(line)
            if m:
                if current_variant is not None and parts:
                    sequences[current_variant] = ''.join(parts).upper()
                current_variant = m.group(1)
                parts = []
                if 'Prototype相同' in line or '与Prototype相同' in line:
                    same_as_proto.append(current_variant)
                    sequences[current_variant] = None
                else:
                    sequences.setdefault(current_variant, '')
                continue
            if SEQUENCE_LINE_PATTERN.match(line):
                seq = re.sub(r'^\s*\d+\s*', '', line)
                seq = re.sub(r'\s+', '', seq)
                parts.append(seq.upper())
    if current_variant is not None and parts:
        sequences[current_variant] = ''.join(parts).upper()

    # Backfill variants marked as Prototype-same
    proto_seq = None
    for key, seq in sequences.items():
        if key and key.lower().startswith('prototype') and seq:
            proto_seq = seq
            break
    if proto_seq:
        for key in same_as_proto:
            sequences[key] = proto_seq

    return {k: v for k, v in sequences.items() if v}


def write_af3_json(job_name: str, sequence: str, out_dir: Path) -> Path:
    payload = {
        "name": job_name,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": sequence}}
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{job_name}_input.json"
    with out_path.open('w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    return out_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data_dir', default='data/按基因分类/按基因分类')
    ap.add_argument('--out_dir', default='af3/input')
    ap.add_argument('--limit', type=int, default=0, help='Optional cap on JSONs written')
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    if not data_dir.exists():
        logger.error("Data directory not found: %s", data_dir)
        return

    total = 0
    for gene_file in sorted(data_dir.glob('*.txt')):
        if gene_file.name.startswith('.'):
            continue
        gene = gene_file.stem
        seqs = parse_gene_file(gene_file)
        if not seqs:
            logger.warning("No sequences parsed: %s", gene_file)
            continue
        for variant, seq in seqs.items():
            job = f"{variant}_{gene}"
            out_path = write_af3_json(job, seq, out_dir)
            total += 1
            logger.info("Wrote %s", out_path)
            if args.limit and total >= args.limit:
                logger.info("Limit reached: %d", args.limit)
                logger.info("Total JSON written: %d", total)
                return
    logger.info("Total JSON written: %d", total)


if __name__ == '__main__':
    main()

