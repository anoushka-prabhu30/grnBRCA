# =============================================================
# Script 08: iDINGO Volcano Plot
# =============================================================
# Generates a volcano plot of differential scores (ds) vs
# -log10(p-value) from iDINGO output.
# Highlights pairs with |ds| > 0.2 and p < 0.05 in red.
# Labels most extreme pairs with |ds| > 4.0 and p < 0.05.
#
# Input:  results/idingo/idingo_unfiltered_result.xlsx
# =============================================================

library(readxl)
library(ggplot2)

df <- read_excel(
  file.path("results", "idingo", "idingo_unfiltered_result.xlsx"),
  sheet = 1
)

df$pair        <- paste(df[[1]], df[[2]], sep = "-")
colnames(df)[5] <- "ds"
colnames(df)[6] <- "pval"

ggplot(df, aes(x = ds, y = -log10(pval))) +
  geom_point(aes(color = abs(ds) > 0.2 & pval < 0.05), size = 2.5) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  geom_text(
    data = subset(df, abs(ds) > 4.0 & pval < 0.05),
    aes(label = pair), vjust = -1.2, size = 3, check_overlap = TRUE
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "blue") +
  labs(
    title = "Volcano Plot of Differential Scores",
    x     = "Differential Score (ds)",
    y     = "-log10(p-value)",
    color = "Significant"
  ) +
  theme_minimal()
