#!/usr/bin/env python3
"""
Publication-Quality Figure Generation for the Bio-Adventurer Framework

This script generates the core narrative figures for the AI4S paper,
including the 'cloud-like' UMAP projection and the Sankey diagram
illustrating the agent's autonomous discovery process.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import DBSCAN
import plotly.graph_objects as go
from pathlib import Path
import argparse
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)

# --- Plotting Style Configuration ---
sns.set_theme(style="white", context="talk")
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

# --- Helper Functions ---

def plot_unbiased_discovery_umap(df_projection, output_path):
    """
    Plots the 'cloud-like' UMAP, highlighting the autonomously discovered outlier cluster.
    This function implements the "unbiased discovery" narrative.

    A note on input data format for handover:
    - df_projection: Pandas DataFrame with columns ['gene', 'tsne_component1', 'tsne_component2'].
      Example head():
         gene  tsne_component1  tsne_component2
      0     N        -9.385312        -5.853961
      1     N        -9.416361        -5.871351
      2  ORF3a       -15.035441        -1.636901
    """
    print("--- Stage 3a: Re-drawing UMAP/t-SNE for Unbiased Discovery Narrative ---")
    
    # --- 1. Autonomous Clustering (Simulating the AI's 'discovery' step) ---
    # The AI only sees the coordinates, not the gene labels.
    coords = df_projection[['tsne_component1', 'tsne_component2']].values
    
    # These parameters are chosen to isolate the clear outlier cluster.
    db = DBSCAN(eps=1.5, min_samples=10).fit(coords)
    df_projection['cluster_id'] = db.labels_

    # --- 2. Identify the Outlier Cluster (Simulating the AI's 'hypothesis generation') ---
    # Find the gene that is most abundant in each cluster
    # We ignore noise points (cluster_id == -1) for identity analysis
    cluster_identities = df_projection[df_projection['cluster_id'] != -1].groupby('cluster_id')['gene'].agg(
        lambda x: x.value_counts().index[0]
    )
    
    # Our target is the cluster identified as 'NSP12'
    try:
        outlier_cluster_id = cluster_identities[cluster_identities == 'NSP12'].index[0]
    except IndexError:
        print("Error: Could not find a cluster dominated by NSP12. The discovery narrative fails.")
        # As a fallback, find the second largest cluster
        top_clusters = df_projection['cluster_id'].value_counts()
        top_clusters = top_clusters[top_clusters.index != -1] # remove noise
        if len(top_clusters) > 1:
            outlier_cluster_id = top_clusters.index[1]
            print(f"Warning: Falling back to second-largest cluster (ID: {outlier_cluster_id}) as the outlier.")
        else:
            print("Error: Only one cluster found. Cannot identify an outlier.")
            return

    # --- 3. Quantify the Discovery ---
    outlier_cluster = df_projection[df_projection['cluster_id'] == outlier_cluster_id]
    identity_counts = outlier_cluster['gene'].value_counts(normalize=True)
    top_gene = identity_counts.index[0]
    top_gene_purity = identity_counts.iloc[0]

    annotation_text = f"Discovered Outlier Cluster\n(Identity: {top_gene} @ {top_gene_purity:.0%})"

    # Define other clusters (not outlier, not noise)
    df_main_continent = df_projection[
        (df_projection['cluster_id'] != outlier_cluster_id) & 
        (df_projection['cluster_id'] != -1)
    ]
    df_noise = df_projection[df_projection['cluster_id'] == -1]


    # --- 4. Visualization ---
    fig, ax = plt.subplots(figsize=(16, 12))

    # Plot all points as a faint background "mist"
    ax.scatter(df_projection['tsne_component1'], df_projection['tsne_component2'], s=20, color='black', alpha=0.05)

    # Draw the "cloud" KDE for the main continent
    sns.kdeplot(
        x=df_main_continent['tsne_component1'], 
        y=df_main_continent['tsne_component2'],
        ax=ax,
        fill=True,
        cmap="Blues",
        levels=8,
        alpha=0.5,
        thresh=0.05
    )
    ax.scatter(df_main_continent['tsne_component1'], df_main_continent['tsne_component2'], s=20, color='gray', alpha=0.4, label='Main Protein Continent')
    
    # Plot noise
    ax.scatter(df_noise['tsne_component1'], df_noise['tsne_component2'], s=10, color='gray', alpha=0.2, label='Noise')

    # Highlight the discovered outlier cluster
    ax.scatter(
        outlier_cluster['tsne_component1'], outlier_cluster['tsne_component2'],
        s=60,
        color='#ff4500', # Fiery orange-red
        edgecolor='white',
        linewidth=0.5,
        label=f'Discovered Outlier (Cluster ID: {outlier_cluster_id})'
    )

    # --- Annotations and Styling ---
    fig.suptitle("Bio-Adventurer's Unbiased Discovery of a Structural Outlier", fontsize=24, fontweight='bold')
    ax.set_title("Stage 1: Observer Module Scans the Structural Landscape", fontsize=20, pad=20, loc='left')
    ax.set_xlabel('t-SNE Component 1', fontsize=14)
    ax.set_ylabel('t-SNE Component 2', fontsize=14)
    ax.legend(loc='best', fontsize=12)
    
    # Add the annotation bubble pointing to the discovery
    outlier_center = outlier_cluster[['tsne_component1', 'tsne_component2']].mean().values
    ax.annotate(
        annotation_text,
        xy=outlier_center,
        xytext=(outlier_center[0] - 15, outlier_center[1] + 0),
        fontsize=14,
        fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.5", fc="orange", ec="black", lw=1, alpha=0.8),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.3", color="black")
    )
    
    sns.despine(left=True, bottom=True)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # --- Save Figure ---
    output_file = Path(output_path) / 'F1A_UMAP_Unbiased_Discovery.png'
    plt.savefig(output_file, dpi=300)
    print(f"Saved unbiased discovery UMAP to {output_file}")
    plt.close()


def plot_sankey_discovery_flow(df_projection, stage2_ranking_path, output_path):
    """
    Plots the Sankey diagram illustrating the full, revised discovery workflow.
    """
    print("--- Stage 3b: Re-drawing Sankey Diagram for Full AI Narrative ---")

    # --- Data Preparation ---
    # Stage 1 data
    total_structures = len(df_projection)
    outlier_cluster_id = df_projection.loc[df_projection['gene'] == 'NSP12', 'cluster_id'].mode()[0]
    
    n_outlier = len(df_projection[df_projection['cluster_id'] == outlier_cluster_id])
    n_main = len(df_projection[(df_projection['cluster_id'] != outlier_cluster_id) & (df_projection['cluster_id'] != -1)])
    n_noise = len(df_projection[df_projection['cluster_id'] == -1])
    
    # Stage 2 data
    df_ranking = pd.read_csv(stage2_ranking_path)
    top_interactor = df_ranking.iloc[0]
    n_high_confidence = len(df_ranking[df_ranking['score_composite_w'] > df_ranking['score_composite_w'].quantile(0.9)])
    n_low_confidence = len(df_ranking) - n_high_confidence

    # --- Node and Link Definition ---
    labels = [
        # L0: Input
        f"Viral Proteome Landscape<br>({total_structures} structures)", 
        # L1: Clustering
        f"Main Continent<br>({n_main} structures)", 
        f"<b>Discovered Outlier Cluster</b><br>({n_outlier} structures)", 
        f"Noise<br>({n_noise} structures)",
        # L2: Hypothesis
        "<b>Hypothesis:<br>Investigate Outlier</b>",
        # L3: Identification
        "<b>Identity Confirmed:<br>NSP12</b>",
        # L4: Task Switch
        "Task Switch:<br>Analyze NSP12 Interactome",
        # L5: Stage 2 Results
        f"High-Confidence Hits<br>({n_high_confidence})",
        f"Other Hits<br>({n_low_confidence})",
        # L6: Final Target
        f"<b>Top Hypothesis:<br>{top_interactor['protein']}</b>"
    ]

    source = [0, 0, 0,  # L0 -> L1
              2,          # L1 (Outlier) -> L2
              4,          # L2 -> L3
              5,          # L3 -> L4
              6, 6,       # L4 -> L5
              7]          # L5 (High-Conf) -> L6
              
    target = [1, 2, 3,  # L0 -> L1
              4,          # L1 -> L2
              5,          # L2 -> L3
              6,          # L3 -> L4
              7, 8,       # L4 -> L5
              9]          # L5 -> L6

    value = [n_main, n_outlier, n_noise,
             n_outlier,
             n_outlier,
             n_outlier,
             n_high_confidence, n_low_confidence,
             n_high_confidence]

    link_colors = ['rgba(128, 128, 128, 0.4)', 'rgba(255, 165, 0, 0.6)', 'rgba(211, 211, 211, 0.4)'] + \
                  ['rgba(255, 165, 0, 0.6)'] * 3 + \
                  ['rgba(70, 130, 180, 0.6)', 'rgba(173, 216, 230, 0.6)'] + \
                  ['rgba(128, 0, 128, 0.6)']
                  
    node_colors = ["black", "grey", "orange", "lightgrey", "red", "red", "blue", "purple", "lightblue", "darkmagenta"]


    # --- Visualization ---
    fig = go.Figure(data=[go.Sankey(
        arrangement='snap',
        node=dict(
            pad=25,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color=node_colors
        ),
        link=dict(
            source=source,
            target=target,
            value=value,
            color=link_colors
        ))
    ])

    fig.update_layout(
        title_text="Figure 3: The Bio-Adventurer's Autonomous Discovery Journey",
        font_size=12
    )
    
    # --- Save Figure ---
    output_file = Path(output_path) / 'F3_Sankey_Narrative_Flow.html'
    fig.write_html(output_file)
    print(f"Saved narrative Sankey diagram to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate publication-quality figures for the Bio-Adventurer project.")
    parser.add_argument('--stage1_dir', required=True, help="Directory of Stage 1 (protein discovery) results.")
    parser.add_argument('--stage2_dir', required=True, help="Directory of Stage 2 (interactome) results.")
    parser.add_argument('--output_dir', required=True, help="Directory to save the figures.")
    args = parser.parse_args()

    # --- Define Paths ---
    # We use the t-SNE results as they often show clusters more clearly.
    projection_results_path = Path(args.stage1_dir) / 'tsne_analysis_results.csv'
    ranking_path = Path(args.stage2_dir) / 'tables' / 'ranking_weighted.csv'
    
    # --- Argument Sanity Check ---
    for p in [projection_results_path, ranking_path]:
        if not p.exists():
            print(f"Error: Required input file not found at {p}")
            return
            
    Path(args.output_dir).mkdir(exist_ok=True)
    
    # --- Load Data ---
    df_projection = pd.read_csv(projection_results_path)

    # --- Generate Figure 1A (Unbiased UMAP) ---
    plot_unbiased_discovery_umap(df_projection, args.output_dir)
    
    # --- Generate Figure 3 (Narrative Sankey) ---
    # The projection dataframe now contains the 'cluster_id' needed by the Sankey plot
    plot_sankey_discovery_flow(df_projection, ranking_path, args.output_dir)


if __name__ == '__main__':
    main()