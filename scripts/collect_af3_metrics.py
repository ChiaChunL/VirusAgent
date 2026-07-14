#!/usr/bin/env python3
import argparse
from pathlib import Path
import json
import re
import pandas as pd


def collect_job(out_dir: Path):
    """Collect a minimal, consistent set of AF3 metrics from a single job folder.

    Preference order for confidence values:
      1) *_summary_confidences.json (AF3 aggregate per job)
      2) *_ranking_scores.csv (if present)
      3) any other *.json that contains numeric keys (fallback; may be noisy)
    """
    metrics = {'job': out_dir.name, 'status': 'unknown'}
    # parse base name
    # strip timestamp prefixes like YYYYMMDD_HHMMSS_
    m = re.match(r"\d+_\d+_(.+)$", out_dir.name)
    if not m:
        m = re.match(r"\d+_(.+)$", out_dir.name)
    base = m.group(1) if m else out_dir.name
    metrics['base'] = base
    # detect files
    files = list(out_dir.glob('**/*'))
    metrics['n_files'] = len(files)
    # Try AF3 summary_confidences first
    scores = {}
    try:
        summaries = list(out_dir.glob('**/*_summary_confidences.json'))
        if summaries:
            obj = json.loads(summaries[0].read_text())
            for k in ['iptm','ptm','ranking_score','fraction_disordered','has_clash']:
                v = obj.get(k, None)
                if isinstance(v, (int,float)):
                    scores[k] = v
    except Exception:
        pass
    # Fallback: parse *_ranking_scores.csv for 'ranking_score'
    if 'ranking_score' not in scores:
        try:
            import csv
            for fp in out_dir.glob('**/*_ranking_scores.csv'):
                with open(fp, 'r') as f:
                    reader = csv.DictReader(f)
                    # take first row's value if exists
                    row = next(reader, None)
                    if row:
                        for k in row:
                            if k.lower().startswith('ranking'):
                                try:
                                    scores['ranking_score'] = float(row[k])
                                except Exception:
                                    pass
                                break
                        break
        except Exception:
            pass
    # Last resort: scan any JSON (may mix fields across files; low priority)
    if not scores:
        for fp in out_dir.glob('**/*.json'):
            try:
                obj = json.loads(fp.read_text())
            except Exception:
                continue
            for k in ['iptm','ptm','ranking_score','confidence','pLDDT','interface_score','ranking_confidence']:
                if isinstance(obj, dict) and k in obj and isinstance(obj[k], (int,float)):
                    scores.setdefault(k, obj[k])
    metrics.update(scores)
    # status heuristic: require summary_confidences or model.cif for 'done'
    has_summary = any(p.name.endswith('_summary_confidences.json') for p in out_dir.glob('**/*'))
    has_model = any(p.name.endswith('.cif') for p in out_dir.glob('**/*'))
    if metrics['n_files'] == 0:
        metrics['status'] = 'empty'
    elif has_summary or has_model:
        metrics['status'] = 'done'
    else:
        metrics['status'] = 'partial'
    return metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--af3_output_root', default='af3/output')
    ap.add_argument('--out_csv', required=True)
    args = ap.parse_args()

    root = Path(args.af3_output_root)
    rows = []
    if root.exists():
        for d in sorted(root.iterdir()):
            if d.is_dir():
                rows.append(collect_job(d))
    df = pd.DataFrame(rows)
    df.to_csv(args.out_csv, index=False)
    print(f'Collected {len(df)} jobs -> {args.out_csv}')

if __name__ == '__main__':
    main()
