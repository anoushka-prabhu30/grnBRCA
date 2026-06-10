# Exploring Gene Regulatory Networks Associated with the Proteasome in Breast Cancer Progression

**Author:** Anoushka A Prabhu (B.Tech Biotechnology)
**Institution:** Bioinformatics Centre (BTIS), Cancer Research Institute, ACTREC — Advanced Centre for Treatment, Research and Education in Cancer, Tata Memorial Centre, Navi Mumbai, India
**Supervisor:** Nikhil Gadewal, Scientific Officer 'F'
**Project Tenure:** January 2 – May 30, 2025

---

## Abstract

Gene regulatory networks (GRNs) are critical for regulating gene expression through interactions between transcription factors and their target genes, driving essential biological processes such as cell growth, differentiation, and responses to environmental signals. Disruptions in GRNs can lead to uncontrolled cell proliferation, evasion of apoptosis, and therapy resistance in cancer.

While extensive research has focused on breast cancer signalling pathways, the regulatory mechanisms surrounding the proteasome remain underexplored. The ubiquitin-proteasome system (UPS) is a crucial component of cellular homeostasis — a multiprotein complex that targets and degrades proteins. Dysregulation of this complex has been linked to carcinogenesis, especially in aggressive subtypes of breast cancer, as tumour suppressors, misfolded proteins, and cell cycle regulators can all be broken down by aberrant proteasomal activity.

This project reconstructs a breast cancer-specific GRN centred on proteasomal components using transcriptome data. By combining three complementary computational approaches — GENIE3, iDINGO, and KiMONo — we identify important regulatory factors and causal relationships that may serve as therapeutic targets or diagnostic biomarkers.

**Keywords:** Multi-Omics, Breast Cancer, Proteasome, Regulators, Network Inference, Systems Biology, Gene Regulatory Network, Prior Knowledge Network

---

## Research Hypothesis

Dysregulation of the ubiquitin-proteasome system (UPS) in breast cancer is driven by specific transcriptional and post-transcriptional regulatory mechanisms that can be systematically identified through integrative gene regulatory network inference using multi-omics data and prior knowledge.

---

## Objectives

1. Utilise a multi-omics framework by integrating gene expression (RNA-seq) and miRNA expression data from TCGA-BRCA samples
2. Employ network inference methods to reconstruct GRNs specific to tumour and normal conditions
3. Incorporate prior biological knowledge (PKN) to improve the specificity and interpretability of the inferred networks
4. Uncover candidate molecular drivers that may serve as potential therapeutic targets or diagnostic biomarkers related to proteasomal dysfunction

---

## Workflow

| Step | Description | Details |
|------|-------------|---------|
| 1 | **Candidate Gene Selection** | 60 genes: 50 proteasomal subunits (PSMA1-7, PSMB1-10, PSMC1-6, PSMD1-14, PSME1-4, PSMG1-4, PSMF1, ADRM1, PAAF1, POMP, USP14) + 10 TFs (NFE2L1, NFE2L2, STAT1, STAT3, FOXO4, TP53, EGR1, NFYA, NFYB, NFYC) |
| 2 | **Prior Knowledge Network (PKN)** | Built via NetworkAnalyst integrating TF-gene, gene-miRNA, TF-miRNA coregulatory and signalling layers → 430 nodes, 696 edges (SIF format) |
| 3 | **Data Acquisition** | TCGA-BRCA via curatedTCGAData (R): RNA-seq (HTSeq-FPKM) + miRNA (RPM), 32 matched tumour-normal pairs (n = 64) |
| 4 | **Data Filtering** | Retain only PKN-relevant genes and miRNAs from the full expression matrices |
| 5 | **Network Inference** | Three parallel approaches: GENIE3, iDINGO, KiMONo |
| 6 | **Key Regulator Identification** | Cross-tool comparison, effect value analysis, hub network ranking |

---

## Repository Structure
```
proteasome-grn-brca/
├── scripts/
│   ├── 01_download_tumor_data.R
│   ├── 02_download_normal_data.R
│   ├── 03_transpose_expression_data.R
│   ├── 04_filter_to_pkn_genes.R
│   ├── 05_genie3_inference.R
│   ├── 06_kimono_tumor_inference.R
│   ├── 07_kimono_normal_inference.R
│   ├── 08_volcano_plot_idingo.R
│   └── 09_kimono_effect_value_plots.R
├── data/
│   ├── pkn/
│   │   ├── PSM_network_analyst.sif
│   │   ├── gene_gene_networkanalyst_kimono_pkn.csv
│   │   └── gene_mirna_networkanalyst_kimono_pkn.csv
│   ├── processed/
│   └── raw/                        (not tracked — regenerate with scripts 01-02)
├── results/
│   ├── genie3/
│   ├── kimono/
│   └── idingo/
├── figures/
├── docs/
│   └── idingo_instructions.md
├── .gitignore
├── LICENSE
└── README.md
```
---

## Pipeline Details

### Script 01 — Download Tumour Data
Uses the `curatedTCGAData` Bioconductor package to download BRCA tumour samples with RNA-seq (`RNASeq2Gene`), miRNA (`miRNASeqGene`), and RPPA assays. Selects 44 samples present across all three assay types to ensure multi-omics consistency. Exports each assay as a CSV file.

### Script 02 — Download Normal Data
Same as Script 01 but extracts adjacent normal tissue samples (TCGA sample code `"10"`). These 32 paired normal samples form the baseline for differential network analysis.

### Script 03 — Transpose Expression Data
Reads the exported CSV files and transposes them from the default format (samples as rows) to the analysis-ready format (genes/miRNAs as rows, samples as columns) using `data.table::transpose`. Produces four transposed CSVs for tumour/normal RNA-seq and miRNA.

### Script 04 — Filter to PKN Genes
Filters the large-scale transposed expression matrices to retain only genes and miRNAs present in the Prior Knowledge Network SIF file exported from NetworkAnalyst. This reduces noise and focuses inference on biologically relevant regulatory elements. Produces `genie3TUMOR.csv`, `genie3NORMAL.csv`, `filtered_mirnas_tumor.csv`, and `filtered_mirnas_normal.csv`.

### Script 05 — GENIE3 Inference
Runs GENIE3 (GEne Network Inference with Ensemble of Trees) using Random Forest regression (1000 trees, K = sqrt, set.seed(123) for reproducibility) on both tumour and normal filtered RNA-seq matrices. Outputs weighted adjacency matrices and ranked link lists. Regulatory edges with weight > 0.065 are used for network visualisation. Also generates regulatory weight heatmaps, weight distribution histograms, and out-degree distributions.

### Script 06 — KiMONo Tumour Inference
Runs KiMONo (Knowledge-guided Multi-Omics Network inference) on tumour RNA-seq and miRNA data jointly, guided by the PKN. KiMONo uses Sparse Group Lasso (SGL) regression to integrate multi-omics data with prior network constraints, producing directed, interpretable networks with effect sizes (regression coefficients) and R-squared model fit values. Uses 2 parallel cores.

### Script 07 — KiMONo Normal Inference
Identical to Script 06 applied to normal tissue samples. Enables direct comparison of tumour vs normal regulatory architectures. Normal models generally exhibited lower R-squared (median approximately 0.16) versus tumour (median approximately 0.25), consistent with greater regulatory heterogeneity in cancer.

### Script 08 — iDINGO Volcano Plot
Generates a volcano plot from iDINGO's differential score output. Points are coloured red when |ds| > 0.2 and p < 0.05. The most extreme pairs (|ds| > 4.0) are labelled. iDINGO itself was run via the Shiny web interface — see docs/idingo_instructions.md.

### Script 09 — KiMONo Effect Value Plots
Produces side-by-side bar plots comparing predictor-wise effect values (regression coefficients) between tumour and normal conditions for selected PSM genes (PSMB8, PSMB9, PSMD9, PSMD10, PSME3) and regulators (NFE2L1, EGR1, NFYA, hsa-let-7b, hsa-let-7c).

---

## Data Sources

| Data Type | Source | Description |
|-----------|--------|-------------|
| RNA-seq | TCGA-BRCA via curatedTCGAData | HTSeq-FPKM gene expression, tumour + normal |
| miRNA-seq | TCGA-BRCA via curatedTCGAData | RPM miRNA expression, tumour + normal |
| PKN | NetworkAnalyst | TF-gene, gene-miRNA, TF-miRNA, signalling interactions |

**Samples:** 32 matched tumour-normal pairs (n = 64 total), selected for availability across all assay types.

Raw data files are not tracked in this repository due to size. Run Scripts 01 and 02 to regenerate.

---

## Prior Knowledge Network (PKN)

The PKN was constructed using NetworkAnalyst (https://www.networkanalyst.ca/) with 50 proteasomal seed genes as input. Four interaction layers were integrated:

1. **TF-Gene Interactions** — from ENCODE, ChEA, and JASPAR
2. **Gene-miRNA Interactions** — from miRTarBase and TargetScan
3. **TF-miRNA Coregulatory Network** — tripartite regulatory structure
4. **Signalling Network** — from KEGG, STRING, and Reactome

**Final PKN stats: 430 nodes, 696 edges** (exported as SIF format, visualised in Cytoscape)

---

## Tools and Methods

| Tool | Version | Method | Reference |
|------|---------|--------|-----------|
| GENIE3 | Bioconductor 1.28.0 | Random Forest ensemble regression | Huynh-Thu et al., 2010 |
| iDINGO | Shiny web app | Partial correlation, chain graph | Class et al., 2018 |
| KiMONo | GitHub cellmapslab | Sparse Group Lasso regression | Ogris et al., 2021 |
| NetworkAnalyst | Web platform | PKN construction | Zhou et al. |
| Cytoscape | Desktop | Network visualisation | |
| curatedTCGAData | Bioconductor | TCGA data access | Ramos et al. |

---

## Key Results

### PKN Summary

| Network | Nodes | Edges |
|---------|-------|-------|
| PKN (before inference) | 430 | 696 |
| Tumour (after KiMONo) | 154 | 581 |
| Normal (after KiMONo) | 172 | 599 |

### Top 10 Regulator-Target Pairs Across All Three Tools

| Rank | GENIE3 (Tumour) | iDINGO (Top Co-Regulator Pairs) | KiMONo (Tumour Predictor-Target) |
|------|----------------|--------------------------------|----------------------------------|
| 1 | PSMB8-PSMB9 | PSMD13 - hsa-mir-103-1 | EGR1 - PSMB3 |
| 2 | EGR1-FOS | PSMD11 - hsa-mir-3677 | PSMD9 - hsa-let-7b |
| 3 | PSMB9-PSMB8 | PSMD4 - hsa-mir-450b | PSMD13 - NFE2L1 |
| 4 | PSME1-PSME2 | PSMD4 - hsa-mir-3166 | PSMD10 - MYC |
| 5 | PSME2-PSME1 | MCM6 - hsa-mir-30c-1 | USF1 - PSMD9 |
| 6 | PSMB4-PSMD4 | POU3F3 - hsa-mir-432 | PSME3 - hsa-let-7b |
| 7 | FOS-EGR1 | POU3F3 - hsa-mir-148a | STAT1 - PSMB9 |
| 8 | ZNF324-ZNF584 | POU3F3 - hsa-mir-3150b | FOS - PSMC2 |
| 9 | PSMA7-ADRM1 | MYC - hsa-mir-181a-2 | PSME3 - NCOA3 |
| 10 | PADI4-NANOG | MYC - hsa-mir-887 | PSMB9 - STAT1 |

### Cross-Tool Consensus Regulators

The following regulators appeared in results from two or more tools, suggesting robust biological relevance:

- **PSMB9** — immunoproteasome subunit; associated with favourable outcomes in triple-negative breast cancer (Geoffroy et al., 2023)
- **EGR1** — tumour-suppressive transcription factor; limits proliferation and promotes apoptosis (Wei et al., 2017)
- **FOS** — restrains metastatic behaviour; role varies by cellular context (Chang et al., 2023)
- **PSMD13** — involved in cellular stemness maintenance and epithelial-mesenchymal transition activation (Holzen et al., 2022)
- **MYC** — major oncogene driving proliferation and genomic instability; poor prognosis across breast cancer subtypes (Qu et al., 2017)

---

## R Package Requirements

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "curatedTCGAData",
  "MultiAssayExperiment",
  "TCGAutils",
  "SummarizedExperiment",
  "GENIE3"
))

# CRAN packages
install.packages(c(
  "data.table", "igraph", "pheatmap", "ggplot2",
  "cowplot", "dplyr", "tidyverse", "DT",
  "oem", "foreach", "doSNOW", "readxl"
))

# KiMONo (from GitHub)
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("cellmapslab/kimono")
```

---

## Future Directions

- Integration of DNA methylation data into the KiMONo framework to incorporate epigenetic regulation
- Expansion of the inference tool suite to include Bayesian, information-theory-based, and correlation-based methods
- Subtype-specific network analysis (Luminal A/B, HER2+, TNBC)
- Wet lab experimental validation of top computational findings

---

## Citation

If you use this work, please cite:

Sundararajan, R., Hegde, S. R., Panda, A. K., Christie, J., Gadewal, N., & Venkatraman, P. (2025). Loss of correlated proteasomal subunit expression selectively promotes the 20S-High state which underlies luminal breast tumorigenicity. *Communications Biology*, 8(1), 55.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
