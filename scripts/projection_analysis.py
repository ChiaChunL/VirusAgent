#!/usr/bin/env python3
"""
AlphaFold3结构预测结果高维投影分析脚本

此脚本从AF3输出结果中提取结构特征，进行高维投影分析，并按变异株和基因分类。
支持多种降维和聚类方法，生成可视化结果。
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import pickle
import warnings
warnings.filterwarnings('ignore')

# 配置中文字体和日志
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AF3ProjectionAnalyzer:
    def __init__(self, af3_results_dir: str, output_dir: str):
        self.af3_results_dir = Path(af3_results_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 创建子目录
        self.embeddings_dir = self.output_dir / "high_dimensional_embeddings"
        self.projections_dir = self.output_dir / "projections"
        self.clusters_dir = self.output_dir / "clusters"
        self.visualizations_dir = self.output_dir / "visualizations"
        
        for dir_path in [self.embeddings_dir, self.projections_dir, 
                        self.clusters_dir, self.visualizations_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
            
        # 变异株颜色映射
        self.variant_colors = {
            'Prototype': '#1f77b4',
            'Alpha': '#ff7f0e', 
            'Beta': '#2ca02c',
            'Gamma': '#d62728',
            'Delta': '#9467bd',
            'Lambda': '#8c564b',
            'Omicron': '#e377c2'
        }
        
    def extract_confidence_features(self, confidence_file: str) -> Optional[np.ndarray]:
        """从AF3置信度文件中提取特征"""
        try:
            with open(confidence_file, 'r') as f:
                confidence_data = json.load(f)
                
            # 提取原子级置信度 (AF3使用atom_plddts字段)
            atom_confidences = confidence_data.get('atom_plddts', [])
            if not atom_confidences:
                # 尝试其他可能的字段名
                atom_confidences = confidence_data.get('confidenceScore', [])
                if not atom_confidences:
                    logger.warning(f"No confidence scores found in {confidence_file}")
                    return None
                
            # 计算统计特征
            atom_conf_array = np.array(atom_confidences)
            features = [
                np.mean(atom_conf_array),
                np.std(atom_conf_array),
                np.median(atom_conf_array),
                np.min(atom_conf_array),
                np.max(atom_conf_array),
                np.percentile(atom_conf_array, 25),
                np.percentile(atom_conf_array, 75),
                len(atom_conf_array)  # 序列长度
            ]
            
            # 添加置信度分布特征
            high_conf = np.sum(atom_conf_array > 90) / len(atom_conf_array)
            medium_conf = np.sum((atom_conf_array > 70) & (atom_conf_array <= 90)) / len(atom_conf_array)
            low_conf = np.sum(atom_conf_array <= 70) / len(atom_conf_array)
            
            features.extend([high_conf, medium_conf, low_conf])
            
            return np.array(features)
            
        except Exception as e:
            logger.error(f"Error extracting features from {confidence_file}: {e}")
            return None
    
    def parse_structure_name(self, structure_name: str) -> Tuple[str, str]:
        """解析结构名称，提取变异株和基因信息"""
        parts = structure_name.split('_')
        if len(parts) >= 2:
            variant = parts[0]
            gene = '_'.join(parts[1:])
            return variant, gene
        return "Unknown", structure_name
    
    def collect_all_features(self) -> Tuple[np.ndarray, List[str], List[str], List[str]]:
        """收集所有AF3结果的特征"""
        features_list = []
        structure_names = []
        variants = []
        genes = []
        
        # 搜索所有confidence文件
        confidence_files = list(self.af3_results_dir.rglob("*_confidences.json"))
        
        logger.info(f"Found {len(confidence_files)} confidence files")
        
        for conf_file in confidence_files:
            # 从文件路径提取结构名称
            structure_name = conf_file.parent.name
            if '_seed-' in structure_name:
                # 移除seed信息，保留主要名称
                structure_name = structure_name.split('_seed-')[0]
                
            # 提取特征
            features = self.extract_confidence_features(conf_file)
            if features is not None:
                variant, gene = self.parse_structure_name(structure_name)
                
                features_list.append(features)
                structure_names.append(structure_name)
                variants.append(variant)
                genes.append(gene)
                
        if not features_list:
            logger.error("No features extracted from any confidence files")
            return None, [], [], []
            
        features_array = np.vstack(features_list)
        logger.info(f"Collected features for {len(features_list)} structures")
        logger.info(f"Feature dimensions: {features_array.shape}")
        
        return features_array, structure_names, variants, genes
    
    def perform_dimensionality_reduction(self, features: np.ndarray, method: str = 'pca') -> np.ndarray:
        """执行降维分析"""
        # 标准化特征
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(features)
        
        logger.info(f"Performing {method.upper()} dimensionality reduction")
        
        if method.lower() == 'pca':
            reducer = PCA(n_components=2)
            projection = reducer.fit_transform(features_scaled)
            
            # 保存解释方差比
            explained_variance = reducer.explained_variance_ratio_
            logger.info(f"PCA explained variance: {explained_variance[0]:.3f}, {explained_variance[1]:.3f}")
            
        elif method.lower() == 'tsne':
            reducer = TSNE(n_components=2, random_state=42, perplexity=min(30, len(features)-1))
            projection = reducer.fit_transform(features_scaled)
            
        elif method.lower() == 'umap':
            try:
                import umap
                reducer = umap.UMAP(n_components=2, random_state=42)
                projection = reducer.fit_transform(features_scaled)
            except ImportError:
                logger.warning("UMAP not available, using PCA instead")
                reducer = PCA(n_components=2)
                projection = reducer.fit_transform(features_scaled)
                
        else:
            logger.warning(f"Unknown method {method}, using PCA")
            reducer = PCA(n_components=2)
            projection = reducer.fit_transform(features_scaled)
            
        # 保存降维结果
        projection_file = self.projections_dir / f"{method}_projection.npy"
        np.save(projection_file, projection)
        
        # 保存降维器
        reducer_file = self.projections_dir / f"{method}_reducer.pkl"
        with open(reducer_file, 'wb') as f:
            pickle.dump((reducer, scaler), f)
            
        return projection
    
    def perform_clustering(self, projection: np.ndarray, variants: List[str], genes: List[str]) -> Dict[str, np.ndarray]:
        """执行聚类分析"""
        clustering_results = {}
        
        # K-means聚类 (基于变异株数量)
        n_variants = len(set(variants))
        kmeans = KMeans(n_clusters=n_variants, random_state=42)
        kmeans_labels = kmeans.fit_predict(projection)
        clustering_results['kmeans'] = kmeans_labels
        
        # DBSCAN聚类
        dbscan = DBSCAN(eps=0.5, min_samples=3)
        dbscan_labels = dbscan.fit_predict(projection)
        clustering_results['dbscan'] = dbscan_labels
        
        # 计算聚类质量
        if len(set(kmeans_labels)) > 1:
            kmeans_score = silhouette_score(projection, kmeans_labels)
            logger.info(f"K-means silhouette score: {kmeans_score:.3f}")
            
        if len(set(dbscan_labels)) > 1:
            dbscan_score = silhouette_score(projection, dbscan_labels)
            logger.info(f"DBSCAN silhouette score: {dbscan_score:.3f}")
            
        # 保存聚类结果
        for method, labels in clustering_results.items():
            cluster_file = self.clusters_dir / f"{method}_clusters.npy"
            np.save(cluster_file, labels)
            
        return clustering_results
    
    def create_visualizations(self, projection: np.ndarray, structure_names: List[str], 
                            variants: List[str], genes: List[str], 
                            clustering_results: Dict[str, np.ndarray], method: str):
        """创建可视化图表"""
        
        # 按变异株着色的散点图
        plt.figure(figsize=(12, 8))
        for variant in set(variants):
            mask = [v == variant for v in variants]
            color = self.variant_colors.get(variant, '#777777')
            plt.scatter(projection[mask, 0], projection[mask, 1], 
                       c=color, label=variant, alpha=0.7, s=60)
        
        plt.xlabel(f'{method.upper()} Component 1')
        plt.ylabel(f'{method.upper()} Component 2')
        plt.title(f'COVID-19 Variants - {method.upper()} Projection')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(self.visualizations_dir / f'{method}_variants.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 按基因着色的散点图
        plt.figure(figsize=(14, 10))
        genes_unique = list(set(genes))
        colors = plt.cm.tab20(np.linspace(0, 1, len(genes_unique)))
        
        for i, gene in enumerate(genes_unique):
            mask = [g == gene for g in genes]
            plt.scatter(projection[mask, 0], projection[mask, 1], 
                       c=[colors[i]], label=gene, alpha=0.7, s=60)
        
        plt.xlabel(f'{method.upper()} Component 1')
        plt.ylabel(f'{method.upper()} Component 2')
        plt.title(f'COVID-19 Genes - {method.upper()} Projection')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
        plt.tight_layout()
        plt.savefig(self.visualizations_dir / f'{method}_genes.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 聚类结果可视化
        for cluster_method, labels in clustering_results.items():
            plt.figure(figsize=(10, 8))
            scatter = plt.scatter(projection[:, 0], projection[:, 1], c=labels, 
                                cmap='viridis', alpha=0.7, s=60)
            plt.colorbar(scatter)
            plt.xlabel(f'{method.upper()} Component 1')
            plt.ylabel(f'{method.upper()} Component 2')
            plt.title(f'{cluster_method.upper()} Clustering on {method.upper()} Projection')
            plt.tight_layout()
            plt.savefig(self.visualizations_dir / f'{method}_{cluster_method}_clustering.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def save_results_dataframe(self, projection: np.ndarray, structure_names: List[str],
                             variants: List[str], genes: List[str],
                             clustering_results: Dict[str, np.ndarray], method: str):
        """保存结果到DataFrame"""
        results_df = pd.DataFrame({
            'structure_name': structure_names,
            'variant': variants,
            'gene': genes,
            f'{method}_component1': projection[:, 0],
            f'{method}_component2': projection[:, 1]
        })
        
        # 添加聚类标签
        for cluster_method, labels in clustering_results.items():
            results_df[f'{cluster_method}_cluster'] = labels
            
        # 保存到CSV
        results_file = self.output_dir / f'{method}_analysis_results.csv'
        results_df.to_csv(results_file, index=False)
        logger.info(f"Results saved to {results_file}")
        
        return results_df
    
    def generate_summary_report(self, results_df: pd.DataFrame, method: str):
        """生成分析总结报告"""
        report_file = self.output_dir / f'{method}_analysis_report.txt'
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(f"AlphaFold3 结构预测结果 {method.upper()} 投影分析报告\n")
            f.write("="*60 + "\n\n")
            
            f.write(f"总结构数量: {len(results_df)}\n")
            f.write(f"变异株数量: {results_df['variant'].nunique()}\n")
            f.write(f"基因数量: {results_df['gene'].nunique()}\n\n")
            
            f.write("变异株分布:\n")
            variant_counts = results_df['variant'].value_counts()
            for variant, count in variant_counts.items():
                f.write(f"  {variant}: {count}\n")
            f.write("\n")
            
            f.write("基因分布:\n")
            gene_counts = results_df['gene'].value_counts()
            for gene, count in gene_counts.items():
                f.write(f"  {gene}: {count}\n")
            f.write("\n")
            
            # 聚类分析
            if 'kmeans_cluster' in results_df.columns:
                f.write("K-means 聚类结果:\n")
                for cluster in sorted(results_df['kmeans_cluster'].unique()):
                    cluster_variants = results_df[results_df['kmeans_cluster'] == cluster]['variant'].value_counts()
                    f.write(f"  聚类 {cluster}: {dict(cluster_variants)}\n")
                f.write("\n")
                
        logger.info(f"Analysis report saved to {report_file}")
    
    def run_full_analysis(self, methods: List[str] = ['pca', 'tsne']):
        """运行完整的投影分析流程"""
        logger.info("Starting AF3 projection analysis...")
        
        # 收集特征
        features, structure_names, variants, genes = self.collect_all_features()
        if features is None:
            logger.error("Failed to collect features")
            return
            
        # 保存原始特征
        features_file = self.embeddings_dir / "raw_features.npy"
        np.save(features_file, features)
        
        metadata_file = self.embeddings_dir / "metadata.json"
        metadata = {
            'structure_names': structure_names,
            'variants': variants,
            'genes': genes,
            'feature_names': [
                'mean_confidence', 'std_confidence', 'median_confidence',
                'min_confidence', 'max_confidence', 'q25_confidence', 'q75_confidence',
                'sequence_length', 'high_conf_ratio', 'medium_conf_ratio', 'low_conf_ratio'
            ]
        }
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
            
        # 对每种降维方法进行分析
        for method in methods:
            logger.info(f"Analyzing with {method.upper()}")
            
            # 降维
            projection = self.perform_dimensionality_reduction(features, method)
            
            # 聚类
            clustering_results = self.perform_clustering(projection, variants, genes)
            
            # 可视化
            self.create_visualizations(projection, structure_names, variants, genes, 
                                     clustering_results, method)
            
            # 保存结果
            results_df = self.save_results_dataframe(projection, structure_names, variants, 
                                                   genes, clustering_results, method)
            
            # 生成报告
            self.generate_summary_report(results_df, method)
            
        logger.info("AF3 projection analysis completed!")

def main():
    """主函数"""
    import argparse
    parser = argparse.ArgumentParser(description="AlphaFold3 结构预测结果高维投影分析脚本")
    parser.add_argument('--af3_results_dir', required=True, help="AF3 aresults directory")
    parser.add_argument('--output_dir', required=True, help="Directory to save analysis results")
    args = parser.parse_args()

    # 检查结果目录是否存在
    if not os.path.exists(args.af3_results_dir):
        logger.error(f"AF3 results directory not found: {args.af3_results_dir}")
        return

    # 创建分析器并运行
    analyzer = AF3ProjectionAnalyzer(args.af3_results_dir, args.output_dir)
    analyzer.run_full_analysis(['pca', 'tsne'])

if __name__ == "__main__":
    main()