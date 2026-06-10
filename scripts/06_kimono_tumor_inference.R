# =============================================================
# Script 06: KiMONo Network Inference — Tumor Condition
# =============================================================
# Runs KiMONo on tumor RNA-seq + miRNA data, guided by PKN.
#
# Input:  data/processed/genie3TUMOR.csv
#         data/processed/filtered_mirnas_tumor.csv
#         data/pkn/gene_gene_networkanalyst_kimono_pkn.csv
#         data/pkn/gene_mirna_networkanalyst_kimono_pkn.csv
# Output: results/kimono/networknew_tumor.csv
# =============================================================

library(kimono)
library(igraph)
library(data.table)
library(oem)
library(foreach)
library(doSNOW)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(DT)

if (!dir.exists(file.path("results", "kimono")))
  dir.create(file.path("results", "kimono"), recursive = TRUE)

transcriptome <- fread(file.path("data", "processed", "genie3TUMOR.csv"))
head(transcriptome)

mirna <- fread(file.path("data", "processed", "filtered_mirnas_tumor.csv"))
head(mirna)

input_data <- list('gene' = transcriptome, 'mirna' = mirna)

gene_gene  <- fread(file.path("data", "pkn", "gene_gene_networkanalyst_kimono_pkn.csv"))
gene_mirna <- fread(file.path("data", "pkn", "gene_mirna_networkanalyst_kimono_pkn.csv"))

prior_network <- create_prior_network(rbind(gene_mirna, gene_gene))

vertex <- do.call(rbind, strsplit(V(prior_network)$name, split = '___'))
prior_network %>% plot(
  edge.curved        = 0,
  main               = 'Prior Network',
  vertex.color       = c("steel blue", "orange")[vertex[, 1] %>% as.factor() %>% as.numeric()],
  vertex.frame.color = "white",
  vertex.label       = vertex[, 2],
  vertex.label.color = 'black',
  vertex.label.cex   = .7,
  layout             = layout_randomly,
  rescale            = F
)
legend(x = -1.5, y = -1.1, c("Genes", "mirna"), pch = 21,
       col = "#777777", pt.bg = c("steel blue", "orange"),
       pt.cex = 2, cex = .8, bty = "n", ncol = 1)

network <- kimono(input_data, prior_network, core = 2, infer_missing_prior = TRUE)

to_igraph(network) %>% plot_kimono(title = 'KiMONo Network (directed)')
to_igraph(network, directed = F) %>% plot_kimono(title = 'KiMONo Network (undirected)')
DT::datatable(head(network), class = 'cell-border stripe')

gg_all <- network[predictor == '(Intercept)', ] %>%
  ggplot(aes(y = r_squared)) + geom_boxplot()

gg_grouped <- network[predictor == '(Intercept)', ] %>%
  ggplot(aes(y = r_squared, x = target_layer, fill = target_layer)) +
  geom_boxplot() +
  scale_fill_manual(values = c("steel blue", "#842F39", "orange"))

plot_grid(gg_all, gg_grouped, rel_widths = c(1, 2))

nnodes <- c(network$target, network$predictor) %>% unique %>% length
nedges <- dim(network)[1]
cat('Number of Nodes: ', nnodes)
cat('Number of Edges: ', nedges)

write.csv(network, file = file.path("results", "kimono", "networknew_tumor.csv"))
message("Saved: results/kimono/networknew_tumor.csv")
