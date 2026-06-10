# =============================================================
# Script 02: Download TCGA-BRCA Normal (Adjacent) Samples
# =============================================================
# Downloads adjacent normal tissue RNA-seq and miRNA data.
# TCGA sample code "10" = solid tissue normal.
#
# Output: data/raw/BRCA_MultiAssayExperiment_Normal/ (CSV files)
# =============================================================

library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)

brca_SE <- curatedTCGAData(
  diseaseCode = "BRCA",
  assays      = c("miRNASeqGene", "RNASeq2Gene", "RPPAArray"),
  version     = "2.0.1",
  dry.run     = FALSE
)

brca_normal <- TCGAsplitAssays(brca_SE, "10")

sampleTables(brca_normal)

export_dir <- file.path("data", "raw", "BRCA_MultiAssayExperiment_Normal")
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

exportClass(brca_normal, dir = export_dir, fmt = "csv", ext = ".csv")

print(paste("Normal sample files exported to:", export_dir))
list.files(export_dir)
