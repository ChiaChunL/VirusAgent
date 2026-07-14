#!/usr/bin/env python3
"""
Narrative Deep Dive Plots for the Bio-Adventurer Framework

This script generates supplementary figures that provide deeper explanations
for the AI's discovery process, including the PCA Loading Plot for Stage 1
and the Feature Importance Bar Chart for Stage 2.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import json
from pathlib import Path
import argparse
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)

# --- Plotting Style Configuration ---
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

def plot_pca_loadings(reducer_path, metadata_path, output_path):
    """
    Visualizes the PCA loadings to explain which features drive the separation.
    """
    print("--- Generating PCA Loading Plot (Explaining Stage 1 'Why') ---")
    
    # --- Load Data ---
    with open(reducer_path, 'rb') as f:
        reducer, scaler = pickle.load(f)
    
    with open(metadata_path, 'r') as f:
        metadata = json.load(f)
        
    feature_names = metadata.get('feature_names', [f'F{i+1}' for i in range(reducer.components_.shape[1])])
    
    # --- Create Plot ---
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Plot the loadings as vectors
    for i, feature in enumerate(feature_names):
        ax.arrow(0, 0, reducer.components_[0, i], reducer.components_[1, i],
                 head_width=0.03, head_length=0.05, fc='red', ec='red', alpha=0.8)
        ax.text(reducer.components_[0, i] * 1.2, reducer.components_[1, i] * 1.2,
                feature, color='black', ha='center', va='center', fontsize=12)

    # --- Styling ---
    ax.set_xlabel('Principal Component 1', fontsize=14)
    ax.set_ylabel('Principal Component 2', fontsize=14)
    ax.set_title('PCA Loading Plot: Contribution of Structural Features', fontsize=18, pad=20)
    
    # Set limits and grid
    lim = np.max(np.abs(reducer.components_)) * 1.5
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axhline(0, color='grey', lw=1, linestyle='--')
    ax.axvline(0, color='grey', lw=1, linestyle='--')
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    
    # --- Save Figure ---
    output_file = Path(output_path) / 'F1B_PCA_Loadings.png'
    plt.savefig(output_file, dpi=300)
    print(f"Saved PCA Loading Plot to {output_file}")
    plt.close()

def plot_feature_weights(weights_path, output_path):
    """
    Visualizes the feature weights from Stage 2 to show what drives the ranking.
    """
    print("--- Generating Feature Weight Bar Chart (Explaining Stage 2 'Why') ---")
    
    # --- Load Data ---
    with open(weights_path, 'r') as f:
        data = json.load(f)
    
    weights = data.get('feature_weights', {})
    if not weights:
        print(f"Warning: 'feature_weights' key not found or empty in {weights_path}. Skipping plot.")
        return

    df_weights = pd.DataFrame(list(weights.items()), columns=['Feature', 'Weight']).sort_values('Weight', ascending=False)
    
    # --- Create Plot ---
    fig, ax = plt.subplots(figsize=(12, 8))
    
    sns.barplot(x='Weight', y='Feature', data=df_weights, palette='viridis', ax=ax)
    
    # --- Styling ---
    ax.set_xlabel('Feature Importance (Entropy Weight)', fontsize=14)
    ax.set_ylabel('Feature', fontsize=14)
    ax.set_title('Stage 2: Objective Feature Weights for Interactor Ranking', fontsize=18, pad=20)
    
    # Add labels to bars
    for i, v in enumerate(df_weights['Weight']):
        ax.text(v + 0.005, i, f"{v:.3f}", color='black', va='center')
        
    plt.tight_layout()
    
    # --- Save Figure ---
    output_file = Path(output_path) / 'F2B_Feature_Weights.png'
    plt.savefig(output_file, dpi=300)
    print(f"Saved Feature Weights chart to {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate deep-dive narrative figures for the Bio-Adventurer project.")
    parser.add_argument('--stage1_dir', required=True, help="Directory of Stage 1 (protein discovery) results.")
    parser.add_argument('--stage2_dir', required=True, help="Directory of Stage 2 (interactome) results.")
    parser.add_argument('--output_dir', required=True, help="Directory to save the figures.")
    args = parser.parse_args()

    # --- Define Paths ---
    pca_reducer_path = Path(args.stage1_dir) / 'projections' / 'pca_reducer.pkl'
    metadata_path = Path(args.stage1_dir) / 'high_dimensional_embeddings' / 'metadata.json'
    feature_weights_path = Path(args.stage2_dir) / 'tables' / 'feature_weights.json'
    
    # --- Argument Sanity Check ---
    for p in [pca_reducer_path, metadata_path, feature_weights_path]:
        if not p.exists():
            print(f"Error: Required input file not found at {p}")
            return
            
    Path(args.output_dir).mkdir(exist_ok=True)

    # --- Generate Plots ---
    plot_pca_loadings(pca_reducer_path, metadata_path, args.output_dir)
    plot_feature_weights(feature_weights_path, args.output_dir)

if __name__ == '__main__':
    main()
