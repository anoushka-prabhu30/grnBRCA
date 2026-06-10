# =============================================================
# Script 05: GENIE3 Gene Regulatory Network Inference
# =============================================================
# Runs GENIE3 (Random Forest GRN inference) on filtered RNA-seq
# for both tumor and normal conditions.
#
# Input:  data/processed/genie3NORMAL.csv
#         data/processed/genie3TUMOR.csv
# Output: results/genie3/genie3normallinklist.csv
#         results/genie3/genie3normalmatrix.csv
#         results/genie3/genie3tumorlinklist.csv
#         results/genie3/genie3tumormatrix.csv
# =============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("GENIE3", quietly = TRUE))
  BiocManager::install("GENIE3")

library(GENIE3)
library(igraph)
library(pheatmap)

if (!dir.exists(file.path("results", "genie3")))
  dir.create(file.path("results", "genie3"), recursive = TRUE)

run_genie3 <- function(csv_path, label) {
  message("\n--- GENIE3: ", label, " ---")
  expr_data <- read.csv(csv_path, row.names = 1)
  exprMatr  <- as.matrix(expr_data)
  message("Matrix: ", nrow(exprMatr), " genes x ", ncol(exprMatr), " samples")

  set.seed(123)
  weightMat <- GENIE3(exprMatr)
  linkList  <- getLinkList(weightMat)

  out <- file.path("results", "genie3", paste0("genie3", tolower(label)))
  write.csv(linkList,  paste0(out, "linklist.csv"), row.names = FALSE)
  write.csv(weightMat, paste0(out, "matrix.csv"),   row.names = TRUE)
  message("Saved: ", out, "linklist.csv and matrix.csv")

  top_links  <- linkList[linkList$weight > 0.065, ]
  g          <- graph_from_data_frame(top_links, directed = TRUE)
  plot(g, edge.arrow.size = 0.3, vertex.size = 10, vertex.label.cex = 0.8,
       main = paste("Regulatory", label, "Gene Network (weight > 0.065)"))

  pheatmap(weightMat,
           color        = colorRampPalette(c("white", "red"))(100),
           main         = paste("Regulatory Weight Heatmap -", label),
           cluster_rows = FALSE, cluster_cols = FALSE,
           fontsize_row = 0.9, fontsize_col = 0.9)

  hist(linkList$weight, breaks = 130,
       main = paste("Distribution of Regulatory Weights -", label),
       xlab = "Weight")

  out_degree <- degree(g, mode = "out")
  hist(out_degree, breaks = 30,
       main = paste("Distribution of Regulatory Out-Degree -", label),
       xlab = "Number of targets per regulator")

  invisible(list(weightMat = weightMat, linkList = linkList, graph = g))
}

normal_results <- run_genie3(file.path("data", "processed", "genie3NORMAL.csv"), "Normal")
tumor_results  <- run_genie3(file.path("data", "processed", "genie3TUMOR.csv"),  "Tumor")
