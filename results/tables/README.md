# Data Tables Documentation

This document provides a summary of the key data tables used for the analysis and visualization of NSP12 interactors. These files represent the original datasets used to generate the reference figures.

---

### `prey_scores.csv`

*   **内容简要 (Summary):** 这是从AP-MS质谱原始数据直接计算出的最详细的得分表。每一行代表一个在某次实验中检测到的“猎物”蛋白，包含了其在实验维度上的多种统计分数。
*   **方法与信息 (Methodology & Information):** 由AP-MS分析流程（如 `03_compute_scores.py`）生成。它通过与阴性对照组比较，计算每个蛋白的富集度（`log2FC_vs_ctrl`）、特异性（`specificity`）、Z-score（`robust_compZ`）和可信度（`WD_score`）等，反映了每个潜在互作的原始实验证据强度。
*   **表头 (Header):**
    ```
    Protein IDs,Protein names,Fasta headers,intensity,strain,bait,is_control,source_file,log2_intensity,ctrl_detect_count,ctrl_log2_median,ctrl_MAD,ctrl_detect_ratio,z_ctrl,log2FC_vs_ctrl,compZ_within_strain,delta_vs_other_baits,robust_compZ,n_detect_baits,n_baits,detect_freq,specificity,WD_score,R,Abundance,MiST_like,promiscuity_score,is_high_conf,is_mid_conf_v2
    ```

---

### `features_matrix.csv`

*   **内容简要 (Summary):** 这是一个特征矩阵，每一行代表一个独特的蛋白质（Bait），列则代表了从AP-MS实验中聚合而来的各种特征。
*   **方法与信息 (Methodology & Information):** 由特征聚合脚本（如 `04_aggregate_features.py`）生成。它将 `prey_scores.csv` 中同一个蛋白质在多次实验中的表现进行汇总（如计算平均值、最大值等），形成一个固定的特征向量，用于后续的机器学习排序。
*   **表头 (Header):**
    ```
    strain,bait,n_prey,n_high,mean_z,max_z,mean_WD,mean_MiST,max_MiST,deg_centrality,eig_centrality,protein,mean_confidence,std_confidence,median_confidence,min_confidence,max_confidence,q25_confidence,q75_confidence,sequence_length,high_conf_ratio,medium_conf_ratio,low_conf_ratio
    ```

---

### `ranking.csv`

*   **内容简要 (Summary):** 一个早期的、包含多种不同排序算法结果的排名表。
*   **方法与信息 (Methodology & Information):** 由无偏排序脚本（如 `05_rank_unbiased.py`）生成。它应用了多种排序算法（如PCA、iForest、Topsis）并给出一个初步的综合得分（`score_composite`），反映了在整合结构信息之前的基线排名。
*   **表头 (Header):**
    ```
    strain,bait,n_prey,n_high,mean_z,max_z,mean_WD,mean_MiST,max_MiST,deg_centrality,eig_centrality,protein,mean_confidence,std_confidence,median_confidence,min_confidence,max_confidence,q25_confidence,q75_confidence,sequence_length,high_conf_ratio,medium_conf_ratio,low_conf_ratio,score_pca1,score_iforest,score_rankagg,score_topsis,score_composite,rank_global,rank_within_strain
    ```

---

### `nsp12_interactors_struct_ifaces.csv`

*   **内容简要 (Summary):** 包含了NSP12与其潜在互作蛋白的AlphaFold3复合体结构预测的关键指标。这是连接实验数据和结构数据的核心摘要表。
*   **方法与信息 (Methodology & Information):** 由AF3结果解析脚本（如 `21_collect_af3_metrics.py`）生成。每一行代表一个NSP12-猎物蛋白复合体，列是各种结构指标，如接口可信度（`struct_iptm`）、接口面积（`iface_area`）、接触对（`iface_contact_pairs`）等。
*   **表头 (Header):**
    ```
    Protein IDs,Protein names,...,struct_iptm,struct_ptm,S_struct,iface_contact_pairs,iface_area,iface_hbonds,iface_salt_bridges,iface_hydroph_contacts
    ```

---

### `nsp12_interactors_struct_ranked.csv`

*   **内容简要 (Summary):** 结合了AP-MS实验证据和AF3结构证据后的主要“无偏见”排名列表。
*   **方法与信息 (Methodology & Information):** 由结构整合排序脚本（如 `24_rank_nsp12_interactors_with_structure.py`）生成。它计算了一个综合了实验和结构证据的最终分数（`score_struct`），并据此进行排序（`priority_rank`），这是项目发现的核心结果。
*   **表头 (Header):**
    ```
    Protein IDs,Protein names,...,S_struct,iface_contact_pairs,...,score_struct,priority_rank
    ```

---

### `nsp12_interactors_struct_ranked_calibrated.csv`

*   **内容简要 (Summary):** 在无偏见排名的基础上，进行“最小先验知识校准”后的最终排名。
*   **方法与信息 (Methodology & Information):** 由校准脚本（如 `27_prior_calibrated_rerank.py`）生成。该步骤旨在将一些已知的、与病毒功能相关的蛋白（如果它们在无偏排名中略低）的权重进行微小提升，得到一个校准后的分数（`score_struct_cal`）和排名（`priority_rank_cal`）。
*   **表头 (Header):**
    ```
    Protein IDs,Protein names,...,score_struct,priority_rank,score_struct_cal,priority_rank_cal
    ```

---

### `pca_analysis_results.csv`

*   **内容简要 (Summary):** 对各个病毒蛋白（非宿主蛋白）的AF3单体结构特征进行主成分分析（PCA）后的结果。
*   **方法与信息 (Methodology & Information):** 由Stage 1的蛋白质发现脚本生成。用于论证为何选择NSP12作为研究对象——因为它在结构特征空间中是一个“异常点”。每一行是一个病毒蛋白，列是其在前两个主成分上的坐标。
*   **表头 (Header):**
    ```
    structure_name,variant,gene,pca_component1,pca_component2,kmeans_cluster,dbscan_cluster
    ```

---

### `nsp12_interface_hotspots.csv`

*   **内容简要 (Summary):** 识别出的NSP12蛋白上与其它蛋白相互作用的“热点”残基。
*   **方法与信息 (Methodology & Information):** 通过分析所有AF3预测的复合体结构，统计NSP12上每个氨基酸残基参与形成界面的频率。高频率的残基即为“热点”。
*   **表头 (Header):**
    ```
    resA_idx,contact_count
    ```

---

### `nsp12_hotspots_conservation.csv`

*   **内容简要 (Summary):** 结合了热点信息与跨物种/病毒株的序列保守性信息。
*   **方法与信息 (Methodology & Information):** 将热点残基的接触频率（`contact_count`）与该残基位点的保守性得分（`conservation_score`）进行关联。用于分析关键的相互作用位点是否在进化中是保守的。
*   **表头 (Header):**
    ```
    position,contact_count,conservation_score,shannon_entropy,domain,functional_site
    ```

---

### `interface_pairs_nsp12.csv`

*   **内容简要 (Summary):** 记录了NSP12与互作蛋白之间具体的残基接触对。
*   **方法与信息 (Methodology & Information):** 从AF3复合体结构中提取的原子级别的详细信息。每一行代表一对相互接触的氨基酸残基（一个来自NSP12，一个来自猎物蛋白）及其距离。
*   **表头 (Header):**
    ```
    job,variant,viral_gene,prey_uniprot,resA_idx,resA_name,resB_idx,resB_name,distance
    ```
