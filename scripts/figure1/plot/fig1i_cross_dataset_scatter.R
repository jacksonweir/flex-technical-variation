#!/usr/bin/env Rscript
# ==============================================================================
# fig1i_cross_dataset_scatter.R
# Cross-dataset scatter: BC002 vs BC011 log2FC in VCC vs PBMC (Figure 1I)
# ==============================================================================
#
# Compares the gene-level BC002 vs BC011 log2FC effect between the VCC (Virtual
# Cell Challenge) dataset and the independent 10x PBMC 16-plex dataset. Points
# are colored by whether the gene is significant (p_val_adj < 0.05) in both
# datasets, one dataset only, or neither.
#
# Significance is called by p_val_adj < 0.05 alone (no |log2FC| threshold).
# Genes that meet significance in both datasets are labeled.
#
# INPUT
#   results/DE_BC002_vs_BC011_lane1_wilcox.csv   — from 02_compute_pairwise_de.R
#   results/DE_BC002_vs_BC011_PBMC_wilcox.csv    — from 07_compute_pbmc_de.R
#
# OUTPUT
#   results/figures/fig1i_scatter_BC002vsBC011_VCC_vs_PBMC.pdf
#
# USAGE
#   Rscript scripts/figure1/plot/fig1i_cross_dataset_scatter.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  library(dplyr)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

PVAL_CUTOFF <- 0.05

VCC_CSV  <- file.path(RESULTS_DIR, "DE_BC002_vs_BC011_lane1_wilcox.csv")
PBMC_CSV <- file.path(RESULTS_DIR, "DE_BC002_vs_BC011_PBMC_wilcox.csv")

for (f in c(VCC_CSV, PBMC_CSV)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f,
         "\nRun the corresponding compute script first.")
  }
}

# ==============================================================================
# Load and merge
# ==============================================================================

df_vcc  <- read.csv(VCC_CSV)
df_pbmc <- read.csv(PBMC_CSV)

df_vcc$p_val_adj[is.na(df_vcc$p_val_adj)]   <- 1
df_pbmc$p_val_adj[is.na(df_pbmc$p_val_adj)] <- 1

df <- merge(
  df_vcc[,  c("gene", "avg_log2FC", "p_val_adj")],
  df_pbmc[, c("gene", "avg_log2FC", "p_val_adj")],
  by = "gene", suffixes = c("_VCC", "_PBMC")
) %>%
  mutate(
    sig_VCC  = p_val_adj_VCC  < PVAL_CUTOFF,
    sig_PBMC = p_val_adj_PBMC < PVAL_CUTOFF,
    category = case_when(
      sig_VCC &  sig_PBMC  ~ "Significant in both",
      sig_VCC & !sig_PBMC  ~ "Significant in VCC only",
      !sig_VCC &  sig_PBMC ~ "Significant in PBMC only",
      TRUE                  ~ "Not significant"
    ),
    category = factor(category, levels = c(
      "Significant in both", "Significant in VCC only",
      "Significant in PBMC only", "Not significant"
    ))
  )

message(sprintf("Genes in common: %d", nrow(df)))
print(table(df$category))

# Spearman correlation on genes significant in either dataset
df_sig <- df[df$sig_VCC | df$sig_PBMC, ]
rho <- cor(df_sig$avg_log2FC_VCC, df_sig$avg_log2FC_PBMC, method = "spearman")
message(sprintf("Spearman rho (genes sig in either dataset): %.3f", rho))

# Label genes significant in both
df$label <- ifelse(df$category == "Significant in both", df$gene, NA_character_)

# ==============================================================================
# Plot
# ==============================================================================

color_vals <- c(
  "Significant in both"     = "tomato3",
  "Significant in VCC only" = "steelblue3",
  "Significant in PBMC only"= "darkorange",
  "Not significant"         = "grey80"
)
size_vals  <- c("Significant in both" = 1.2, "Significant in VCC only" = 0.8,
                "Significant in PBMC only" = 0.8, "Not significant" = 0.4)
alpha_vals <- c("Significant in both" = 0.9, "Significant in VCC only" = 0.7,
                "Significant in PBMC only" = 0.7, "Not significant" = 0.3)

p <- ggplot(df[order(df$category, decreasing = TRUE), ],
            aes(x = avg_log2FC_VCC, y = avg_log2FC_PBMC,
                color = category, size = category, alpha = category)) +
  geom_point(shape = 16) +
  geom_text_repel(
    aes(label = label),
    size = 3, color = "black", max.overlaps = 15,
    box.padding = 0.2, point.padding = 0.1,
    segment.size = 0.3, segment.color = "grey50",
    show.legend = FALSE, na.rm = TRUE
  ) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_color_manual(values = color_vals) +
  scale_size_manual(values  = size_vals) +
  scale_alpha_manual(values = alpha_vals) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("Spearman \u03c1 = %.3f", rho),
           size = 4, color = "black") +
  labs(
    x = "Average log2 Fold Change (BC002 / BC011, VCC lane 1)",
    y = "Average log2 Fold Change (BC002 / BC011, PBMC)"
  ) +
  theme_few() +
  theme(
    axis.title        = element_text(size = 17, color = "black"),
    axis.text         = element_text(size = 14, color = "black"),
    legend.position   = "none",
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "fig1i_scatter_BC002vsBC011_VCC_vs_PBMC.pdf"),
       p, width = 7, height = 5, device = cairo_pdf)
message("Saved: fig1i_scatter_BC002vsBC011_VCC_vs_PBMC.pdf")

message("Done.")
print(sessionInfo())
