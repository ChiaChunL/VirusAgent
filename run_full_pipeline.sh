#!/usr/bin/env bash
set -euo pipefail

# This script orchestrates the full two-stage analysis pipeline.
# Stage 1: Unbiased discovery of key viral proteins using PCA/t-SNE on AlphaFold3 structural embeddings.
# Stage 2: Unbiased discovery of host interactors for the key viral protein (NSP12).

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
AF3_RESULTS_DIR="$ROOT_DIR/data/af3_results"
STAGE1_OUTPUT_DIR="$ROOT_DIR/results/stage1_protein_discovery"
STAGE2_OUTPUT_DIR="$ROOT_DIR/results/stage2_interactome"
FIGS_DIR="$ROOT_DIR/results/figs"

# Create the figures directory if it doesn't exist
mkdir -p "$FIGS_DIR"

echo "--- Starting Stage 1: Unbiased Discovery of Key Viral Proteins ---"

python "$ROOT_DIR/scripts/projection_analysis.py" \
    --af3_results_dir "$AF3_RESULTS_DIR" \
    --output_dir "$STAGE1_OUTPUT_DIR"

echo "--- Stage 1 Complete. Analysis saved in $STAGE1_OUTPUT_DIR ---"
echo ""
echo "--- Starting Stage 2: Unbiased Discovery of Host Interactors for NSP12 ---"

bash "$ROOT_DIR/run_stage2_interactome.sh"

echo "--- Stage 2 Complete. ---"
echo ""
echo "--- Stage 3: Publication Figure Generation ---"
# Ensure necessary python packages are installed
pip install -q umap-learn plotly

# Generate main publication figures (UMAP, Sankey)
python3 scripts/publication_figures.py \
    --stage1_dir "$STAGE1_OUTPUT_DIR" \
    --stage2_dir "$STAGE2_OUTPUT_DIR" \
    --output_dir "$FIGS_DIR"
echo "--- Publication figures generated in $FIGS_DIR ---"

# Generate deep-dive explanatory figures (PCA Loadings, Feature Weights)
python3 scripts/narrative_deep_dive_plots.py \
    --stage1_dir "$STAGE1_OUTPUT_DIR" \
    --stage2_dir "$STAGE2_OUTPUT_DIR" \
    --output_dir "$FIGS_DIR"
echo "--- Deep-dive figures generated in $FIGS_DIR ---"

echo ""
echo "Full pipeline finished successfully."
