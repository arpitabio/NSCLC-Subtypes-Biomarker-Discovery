# NSCLC-Subtypes-Biomarker-Discovery
Comparative transcriptomic and network-based analysis of LUAD and LUSC from TCGA to identify subtype-specific and shared biomarkers in non-small cell lung cancer (NSCLC).
# 🧬 NSCLC Biomarker Discovery using TCGA: LUAD and LUSC

This repository contains R scripts for transcriptomic and network-based analysis of non-small cell lung cancer (NSCLC), specifically focusing on the two major subtypes: **lung adenocarcinoma (LUAD)** and **lung squamous cell carcinoma (LUSC)**.

The workflow includes RNA-Seq data acquisition, differential expression analysis, functional enrichment, and biomarker discovery through network topology.

---

## 📌 Objectives

- Retrieve and preprocess TCGA RNA-Seq data for LUAD and LUSC  
- Perform differential gene expression analysis using DESeq2  
- Conduct GO, KEGG, and Disease Ontology enrichment (ORA & GSEA)  
- Visualize PCA, volcano plots, MA plots, dot plots, ridge plots, and enrichment maps  
- Identify hub genes using Cytoscape with centrality measures  
- Apply MCC to prioritize final biomarkers  
- Compare LUAD and LUSC to extract shared molecular signatures

---

## 🧪 Workflow Summary

### 1. **Data Acquisition and Preprocessing**
- Downloaded LUAD and LUSC RNA-Seq data and clinical metadata from TCGA via `TCGAbiolinks`
- Verified sample types and filtered low-expression genes
- Applied log2(count + 1) transformation and performed PCA for quality control

### 2. **Differential Gene Expression Analysis**
- DESeq2 used for tumor vs. normal contrast
- Filtered genes based on padj < 0.05 and |log2FC| > 1
- Annotated gene IDs using `biomaRt` and `org.Hs.eg.db`

### 3. **Functional Enrichment**
- ORA: GO (BP, MF, CC), KEGG, and Disease Ontology terms
- GSEA: Pre-ranked analysis using log2FC with `gseGO` and `gseKEGG`
- Visualized top enriched categories with dot plots, bar plots, ridge plots, and enrichment maps

### 4. **Biomarker Discovery and Network Analysis**
- Upregulated and downregulated DEGs analyzed using Cytoscape (degree, closeness, betweenness)
- MCC (Maximal Clique Centrality) scoring used to identify final biomarkers
- Compared LUAD and LUSC biomarkers to extract shared core genes

### 5. **Validation and Visualization**
- Applied Benjamini–Hochberg correction for multiple testing
- Cross-referenced candidate biomarkers with literature
- Exported figures: volcano plots, MA plots, GSEA plots, network diagrams

---

## 🛠️ Tools & Packages

- **R (v4.4.3)**  
  - `TCGAbiolinks`, `DESeq2`, `clusterProfiler`, `biomaRt`, `org.Hs.eg.db`, `ggplot2`, `dplyr`, etc.  
- **Cytoscape (v3.10.3)**  
  - Plugins: `stringApp (v2.2.0)`, `cytoHubba (v0.1)`

  ## References
- TCGA Data Portal: https://portal.gdc.cancer.gov
- DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
- clusterProfiler: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
- Cytoscape: https://cytoscape.org
  
## 📫 Contact

For questions or suggestions, feel free to reach out via [gangulyarpita561@gmail.com].

---

