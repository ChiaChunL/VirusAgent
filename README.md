# 🦠 VirusAgent

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![AlphaFold3](https://img.shields.io/badge/AlphaFold3-6DB33F.svg)](https://github.com/google-deepmind/alphafold3)
[![Status](https://img.shields.io/badge/manuscript-under%20review-orange.svg)](#citation)

**VirusAgent** is an autonomous scientific discovery pipeline for SARS-CoV-2 host–pathogen interactome analysis. It integrates AP-MS proteomics, structural prediction, and multi-modal evidence calibration through a three-stage **Observe–Investigate–Critique** cognitive cycle.

---

## ⚙️ Installation

**Requirements:** Python 3.9+, pip.

```bash
git clone https://github.com/ChiaChunL/VirusAgent.git
cd VirusAgent
pip install numpy pandas scipy scikit-learn matplotlib adjusttext biopython
```

AlphaFold3 must be installed separately for Stage 2 structural predictions. See the [AlphaFold3 repository](https://github.com/google-deepmind/alphafold3) for setup instructions.

---

## 🚀 Quick Start

**1. Configure local paths:**

```bash
cp configs/config.yaml.example configs/config.yaml
# Edit config.yaml to point to your local data directories
```

**2. Run the full pipeline (Stages 1–3):**

```bash
bash run_full_pipeline.sh
```

**3. Or run each stage individually:**

```bash
# Stage 1 — Observer: AP-MS normalization and viral protein outlier detection
python scripts/projection_analysis.py

# Stage 2 — Investigator: AP-MS scoring and interactor ranking
bash run_stage2_interactome.sh

# Stage 3 — Critic: MORPH calibration and final ranking
python scripts/prior_calibrated_rerank.py
python scripts/storyline_figures.py
```

---

## 🗂️ Repository Layout

```
VirusAgent/
├── configs/
│   └── config.yaml.example                    # Copy to config.yaml and set local paths
├── scripts/
│   ├── standardize_apms.py                    # Stage 1: AP-MS normalization
│   ├── build_control_background.py            # Stage 1: control background
│   ├── compute_scores.py                      # Stage 1: robust Z-score & specificity
│   ├── aggregate_features.py                  # Stage 1: multi-feature matrix
│   ├── rank_unbiased.py                       # Stage 1: unbiased candidate ranking
│   ├── compute_entropy_weights.py             # Stage 1: entropy-based feature weighting
│   ├── rank_weighted.py                       # Stage 1: weighted ranking
│   ├── projection_analysis.py                 # Stage 1: PCA/t-SNE viral protein outlier detection
│   ├── build_nsp12_af3_complex_requests.py    # Stage 2: AlphaFold3 job construction
│   ├── rank_nsp12_interactors_with_structure.py  # Stage 2: structure-informed ranking
│   ├── prior_calibrated_rerank.py             # Stage 3: MORPH calibration
│   └── storyline_figures.py                   # Stage 3: final publication figures
├── results/
│   ├── figs/                                  # Publication figures
│   └── tables/                                # Final ranked interactor tables
├── run_full_pipeline.sh                       # Master pipeline script (Stages 1–3)
└── run_stage2_interactome.sh                  # Stage 2 only
```

---

## 🔬 Pipeline Overview

| Stage | Name | Input | Key Steps | Output |
|-------|------|-------|-----------|--------|
| **1** | Observer | Multi-strain AP-MS data | Normalization → Robust Z-score → PCA/t-SNE | NSP12 flagged as viral bait |
| **2** | Investigator | Top-50 AP-MS candidates | AlphaFold3 → ipTM + interface area scoring | Structural scores for 30 interactors |
| **3** | Critic | Stage 1–2 scores | MORPH fusion → Dirichlet perturbation (500 trials) → ε-margin calibration | Final calibrated ranking |

---

## 📄 Citation

If you use VirusAgent in your research, please cite:

```bibtex
@article{virusagent2026,
  title   = {AI-driven interactome analysis reveals the immune evasion by SARS-CoV-2 RdRp counteracting RIG-I-mediated antiviral activity},
  journal = {Manuscript under review},
  year    = {2026}
}
```

*(This entry will be updated with the final journal, volume, and DOI upon publication.)*

---

## 🙏 Acknowledgements

Structural predictions were performed using [AlphaFold3](https://github.com/google-deepmind/alphafold3) ([Jumper et al., *Nature*, 2024](https://www.nature.com/articles/s41586-024-07487-w)). We thank the AlphaFold team for making their model publicly available.

---

## 📜 License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

---

## ✉️ Contact

For inquiries and collaboration opportunities, feel free to contact us at ChiaChun.Le@gmail.com.
