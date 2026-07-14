#!/usr/bin/env python3
"""
One-shot collector to copy completed AF3 outputs from af3/output/*_<job>
into af3_results/by_variant and af3_results/by_gene.

Criteria for completion: a *_confidences.json exists within the job dir.

Usage:
  python final/scripts/03_collect_af3_results_once.py \
    --af3_output_root af3/output \
    --results_root af3_results
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("collect_af3")


def discover_jobs(af3_output_root: Path) -> list[tuple[str, Path]]:
    jobs: list[tuple[str, Path]] = []
    for ts_dir in sorted(af3_output_root.glob('*_*')):
        if not ts_dir.is_dir():
            continue
        for job_dir in ts_dir.iterdir():
            if not job_dir.is_dir():
                continue
            job = job_dir.name
            jobs.append((job, ts_dir))
    return jobs


def is_completed(job_root: Path, job: str) -> bool:
    # Accept either aggregate or per-sample confs
    if (job_root / job / f"{job}_confidences.json").exists():
        return True
    for p in (job_root / job).glob('seed-*_*/*_confidences.json'):
        if p.exists():
            return True
    return False


def copy_results(job: str, variant: str, gene: str, source_ts_dir: Path, dest_root: Path) -> None:
    by_variant = dest_root / 'by_variant' / f"{variant}_genes" / job
    by_gene = dest_root / 'by_gene' / f"{gene}_variants" / job
    for dest in (by_variant, by_gene):
        if dest.exists():
            shutil.rmtree(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source_ts_dir, dest)
        logger.info("Copied -> %s", dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--af3_output_root', default='af3/output')
    ap.add_argument('--results_root', default='af3_results')
    args = ap.parse_args()

    out_root = Path(args.af3_output_root)
    res_root = Path(args.results_root)
    res_root.mkdir(parents=True, exist_ok=True)

    jobs = discover_jobs(out_root)
    logger.info("Found %d timestamped job roots", len(jobs))

    done = 0
    for job, ts_dir in jobs:
        parts = job.split('_', 1)
        if len(parts) < 2:
            continue
        variant, gene = parts[0], parts[1]
        job_root = ts_dir / job
        if not job_root.exists():
            # Sometimes outputs are directly under timestamp dir
            job_root = ts_dir
        if is_completed(job_root, job):
            copy_results(job, variant, gene, ts_dir, res_root)
            done += 1
    logger.info("Copied %d completed jobs", done)


if __name__ == '__main__':
    main()

