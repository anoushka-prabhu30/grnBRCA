# =============================================================
# Script 03: Transpose RNA-seq and miRNA Expression Data
# =============================================================
# Reads exported CSVs from scripts 01/02.
# Transposes so rows = genes/miRNAs, columns = samples.
#
# NOTE: Adjust filenames to match what curatedTCGAData exported.
# =============================================================

library(data.table)

transpose_and_save <- function(input_path, output_path) {
  dt <- fread(input_path)
  dt$rownames_col <- rownames(dt)
  dt_transposed <- transpose(dt, keep.names = "V1")
  if (!dir.exists(dirname(output_path))) {
    dir.create(dirname(output_path), recursive = TRUE)
  }
  fwrite(dt_transposed, output_path)
  message("Saved: ", output_path)
}

# Normal samples
transpose_and_save(
  input_path  = file.path("data", "raw", "BRCA_MultiAssayExperiment_Normal",
                           "rnaseq_normal.csv"),
  output_path = file.path("data", "processed", "RNASeqnormal_transposed.csv")
)
transpose_and_save(
  input_path  = file.path("data", "raw", "BRCA_MultiAssayExperiment_Normal",
                           "mirna_normal.csv"),
  output_path = file.path("data", "processed", "mirnanormal_transposed.csv")
)

# Tumor samples
transpose_and_save(
  input_path  = file.path("data", "raw", "BRCA_MultiAssayExperiment_Tumor",
                           "rnaseq_tumor.csv"),
  output_path = file.path("data", "processed", "brca_RNASeqfilter_transposed.csv")
)
transpose_and_save(
  input_path  = file.path("data", "raw", "BRCA_MultiAssayExperiment_Tumor",
                           "mirna_tumor.csv"),
  output_path = file.path("data", "processed", "brca_miRNASeqGene_transposed.csv")
)
