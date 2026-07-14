#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
CFG="$ROOT/configs/config.yaml"
# tiny yaml reader via grep/sed (avoid extra deps)
get_cfg(){ local key=$1; grep -E "^\s*${key}:" "$CFG" | sed -E "s/^[^:]+:\s*//"; }
APMS_ROOT=$(get_cfg 'apms_root')
CONTROLS_DIR=$(get_cfg 'controls_dir')
SEQS_DIR=$(get_cfg 'sequences_by_gene')
HIGHDIM_DIR=$(get_cfg 'highdim_dir')
INTER_DIR=$(get_cfg 'intermediate_dir')
TABLES_DIR=$(get_cfg 'tables_dir')
FIGS_DIR=$(get_cfg 'figs_dir')

mkdir -p "$INTER_DIR" "$TABLES_DIR" "$FIGS_DIR"

echo "[1/9] Standardize AP-MS -> apms_long.csv"
python "$ROOT/scripts/standardize_apms.py" \
  --apms_root "$APMS_ROOT" \
  --out_csv "$INTER_DIR/apms_long.csv"

echo "[2/9] Control background stats -> control_stats.csv"
python "$ROOT/scripts/build_control_background.py" \
  --apms_long "$INTER_DIR/apms_long.csv" \
  --out_csv "$INTER_DIR/control_stats.csv" \
  --pseudo_percentile 5

echo "[3/9] Compute robust & specificity scores -> prey_scores.csv"
python "$ROOT/scripts/compute_scores.py" \
  --apms_long "$INTER_DIR/apms_long.csv" \
  --control_stats "$INTER_DIR/control_stats.csv" \
  --out_csv "$TABLES_DIR/prey_scores.csv" \
  --z_thr 3.0 --log2fc_thr 1.0 --pseudo_percentile 5

echo "[4/9] Aggregate features + merge structure -> features_matrix.csv"
python "$ROOT/scripts/aggregate_features.py" \
  --scores_csv "$TABLES_DIR/prey_scores.csv" \
  --highdim_dir "$HIGHDIM_DIR" \
  --out_csv "$TABLES_DIR/features_matrix.csv"

echo "[5/9] Unbiased ranking & evaluation -> ranking.csv"
python "$ROOT/scripts/rank_unbiased.py" \
  --features_csv "$TABLES_DIR/features_matrix.csv" \
  --out_csv "$TABLES_DIR/ranking.csv" \
  --eval_txt "$ROOT/logs/nsp12_eval.txt"

echo "[6/9] Compute objective (entropy) feature weights -> feature_weights.json"
python "$ROOT/scripts/compute_entropy_weights.py" \
  --features_csv "$TABLES_DIR/features_matrix.csv" \
  --out_json "$TABLES_DIR/feature_weights.json"

echo "[7/9] Weighted ranking & evaluation -> ranking_weighted.csv"
python "$ROOT/scripts/rank_weighted.py" \
  --features_csv "$TABLES_DIR/features_matrix.csv" \
  --weights_json "$TABLES_DIR/feature_weights.json" \
  --out_csv "$TABLES_DIR/ranking_weighted.csv" \
  --eval_txt "$ROOT/logs/nsp12_eval_weighted.txt"

echo "[8/9] Visualizations -> figs"
python "$ROOT/scripts/visualize_rankings.py" \
  --ranking_csv "$TABLES_DIR/ranking.csv" \
  --figs_dir "$FIGS_DIR"

echo "[9/9] Advanced Visualizations -> figs"
pip install -q adjustText
python "$ROOT/scripts/advanced_visualize.py" \
  --prey_scores_csv "$TABLES_DIR/prey_scores.csv" \
  --figs_dir "$FIGS_DIR"

echo "Pipeline complete. See $TABLES_DIR and $FIGS_DIR"
