# =============================================================
# Script 04: Filter RNA-seq and miRNA Data to PKN-Relevant Genes
# =============================================================
# Keeps only genes/miRNAs present in the Prior Knowledge Network
# exported from NetworkAnalyst. Applied to tumor and normal.
#
# Input:  data/processed/ (transposed CSVs)
#         data/pkn/PSM_network_analyst.sif
# Output: data/processed/genie3TUMOR.csv
#         data/processed/genie3NORMAL.csv
#         data/processed/filtered_mirnas_tumor.csv
#         data/processed/filtered_mirnas_normal.csv
# =============================================================

library(data.table)

# Load PKN SIF file — 3 columns: source, interaction_type, target
pkn <- fread(file.path("data", "pkn", "PSM_network_analyst.sif"),
             header = FALSE, col.names = c("source", "type", "target"))

pkn_nodes <- unique(c(pkn$source, pkn$target))
message("Total PKN nodes: ", length(pkn_nodes))

filter_to_pkn <- function(input_path, output_path, nodes) {
  data <- fread(input_path)
  filtered <- data[data[[1]] %in% nodes, ]
  if (!dir.exists(dirname(output_path))) {
    dir.create(dirname(output_path), recursive = TRUE)
  }
  fwrite(filtered, output_path)
  message("Filtered: ", output_path, " | Rows kept: ", nrow(filtered))
}

filter_to_pkn(
  file.path("data", "processed", "brca_RNASeqfilter_transposed.csv"),
  file.path("data", "processed", "genie3TUMOR.csv"),
  pkn_nodes
)
filter_to_pkn(
  file.path("data", "processed", "RNASeqnormal_transposed.csv"),
  file.path("data", "processed", "genie3NORMAL.csv"),
  pkn_nodes
)
filter_to_pkn(
  file.path("data", "processed", "brca_miRNASeqGene_transposed.csv"),
  file.path("data", "processed", "filtered_mirnas_tumor.csv"),
  pkn_nodes
)
filter_to_pkn(
  file.path("data", "processed", "mirnanormal_transposed.csv"),
  file.path("data", "processed", "filtered_mirnas_normal.csv"),
  pkn_nodes
)

message("Filtering complete.")
