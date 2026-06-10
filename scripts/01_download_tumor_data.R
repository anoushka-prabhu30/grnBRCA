# =============================================================
# Script 01: Download TCGA-BRCA Tumor Samples
# =============================================================
# Downloads matched tumor RNA-seq and miRNA data from TCGA-BRCA
# using curatedTCGAData. Finds samples present in all 3 assays,
# subsets to 44 common cases, and exports as CSV files.
#
# Output: data/raw/BRCA_MultiAssayExperiment_Tumor/ (CSV files)
# Run once — large download (~1-2 GB).
# =============================================================

library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(SummarizedExperiment)

brca_SE <- curatedTCGAData(
  diseaseCode = "BRCA",
  assays      = c("miRNASeqGene", "RNASeq2Gene", "RPPAArray"),
  version     = "2.0.1",
  dry.run     = FALSE
)

brca_SE <- TCGAprimaryTumors(brca_SE)

random_samples <- sample(
  names(which(table(sampleMap(brca_SE)$primary) == 3)), 44
)
brca_SE <- subsetByColData(brca_SE, random_samples)

export_dir <- file.path("data", "raw", "BRCA_MultiAssayExperiment_Tumor")
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

exportClass(brca_SE, dir = export_dir, fmt = "csv", ext = ".csv")

print(paste("Exported files stored in:", export_dir))
list.files(export_dir)
