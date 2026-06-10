# =============================================================
# Script 09: KiMONo Effect Value Bar Plots
# =============================================================
# Compares predictor-wise effect values (regression coefficients)
# between tumor and normal conditions for selected PSM genes,
# TFs, and miRNAs from KiMONo output.
#
# Input:  results/kimono/networknew_tumor.csv
#         results/kimono/normal_kimono_result_new_pkn.csv
# =============================================================

library(ggplot2)
library(dplyr)
library(data.table)
library(cowplot)

tumor_net  <- fread(file.path("results", "kimono", "networknew_tumor.csv"))
normal_net <- fread(file.path("results", "kimono", "normal_kimono_result_new_pkn.csv"))

plot_effect_values <- function(target_gene, tumor_net, normal_net,
                               exclude_predictors = NULL) {
  t_data <- tumor_net[target == target_gene & predictor != "(Intercept)"]
  n_data <- normal_net[target == target_gene & predictor != "(Intercept)"]

  if (!is.null(exclude_predictors)) {
    t_data <- t_data[!predictor %in% exclude_predictors]
    n_data <- n_data[!predictor %in% exclude_predictors]
  }

  t_data$Type <- "Tumor"
  n_data$Type <- "Normal"
  combined    <- rbind(t_data, n_data)

  if (nrow(combined) == 0) { message("No data for: ", target_gene); return(NULL) }

  ggplot(combined, aes(x = predictor, y = value, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Normal" = "#F08080", "Tumor" = "#20B2AA")) +
    labs(
      title = paste("Effector Values for", target_gene, ": Normal vs Tumor"),
      x = "Predictor", y = "Effector Value"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

psm_genes <- c("PSMB8", "PSMB9", "PSMD9", "PSMD10", "PSME3")
psm_plots <- lapply(psm_genes, plot_effect_values,
                    tumor_net = tumor_net, normal_net = normal_net)
psm_plots <- Filter(Negate(is.null), psm_plots)
if (length(psm_plots) > 0) do.call(plot_grid, c(psm_plots, list(ncol = 3)))

reg_genes <- c("NFE2L1", "EGR1", "NFYA", "hsa-let-7b", "hsa-let-7c")
reg_plots <- lapply(reg_genes, function(g) {
  excl <- if (g == "hsa-let-7b") "PSME3" else NULL
  plot_effect_values(g, tumor_net, normal_net, exclude_predictors = excl)
})
reg_plots <- Filter(Negate(is.null), reg_plots)
if (length(reg_plots) > 0) do.call(plot_grid, c(reg_plots, list(ncol = 3)))
