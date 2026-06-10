# Exploring Gene Regulatory Networks Associated with the Proteasome in Breast Cancer Progression

**Author:** Anoushka A Prabhu
**Institution:** BTIS, CRI, ACTREC, Navi Mumbai
**Supervisor:** Nikhil Gadewal, Scientific Officer F
**Project Tenure:** January 2 - May 30, 2025

## Overview
This project investigates the gene regulatory network (GRN) surrounding the ubiquitin-proteasome system (UPS) in breast cancer using multi-omics data from TCGA-BRCA. Three complementary network inference approaches were applied to identify key transcriptional and post-transcriptional regulators linked to proteasomal dysfunction.

## Pipeline
| Script | Description |
|--------|-------------|
| 01_download_tumor_data.R | Download TCGA-BRCA tumor samples via curatedTCGAData |
| 02_download_normal_data.R | Download TCGA-BRCA normal samples |
| 03_transpose_expression_data.R | Transpose RNA-seq and miRNA matrices |
| 04_filter_to_pkn_genes.R | Filter to PKN-relevant genes and miRNAs |
| 05_genie3_inference.R | GENIE3 Random Forest GRN inference |
| 06_kimono_tumor_inference.R | KiMONo multi-omics inference - tumor |
| 07_kimono_normal_inference.R | KiMONo multi-omics inference - normal |
| 08_volcano_plot_idingo.R | iDINGO volcano plot |
| 09_kimono_effect_value_plots.R | KiMONo effect value comparison plots |

## Tools
- GENIE3 (Bioconductor) - Random Forest network inference
- iDINGO (Shiny web app) - Differential partial correlation network
- KiMONo (R) - Sparse Group Lasso multi-omics inference
- NetworkAnalyst - Prior Knowledge Network construction
- TCGA-BRCA dataset via curatedTCGAData

## Key Findings
Top regulators identified across all 3 tools include EGR1, PSMB8/PSMB9, NFE2L1, MYC, STAT1, and hsa-let-7b as key drivers of proteasomal dysregulation in breast cancer.
