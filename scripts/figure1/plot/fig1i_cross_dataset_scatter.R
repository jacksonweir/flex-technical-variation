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
# x-axis = PBMC log2FC, y-axis = VCC log2FC.
# Dashed lines at ±0.3 log2FC on both axes.
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
# PACKAGES
#   ggplot2, ggrepel, dplyr
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

PVAL_CUTOFF <- 0.05
FC_LINE     <- 0.3

VCC_CSV  <- file.path(RESULTS_DIR, "DE_BC002_vs_BC011_lane1_wilcox.csv")
PBMC_CSV <- file.path(RESULTS_DIR, "DE_BC002_vs_BC011_PBMC_wilcox.csv")

for (f in c(VCC_CSV, PBMC_CSV)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f,
         "\nRun the corresponding compute script first.")
  }
}

# ==============================================================================
# Load and merge (x = PBMC, y = VCC)
# ==============================================================================

df_pbmc <- read.csv(PBMC_CSV)
df_vcc  <- read.csv(VCC_CSV)

if (!("gene" %in% colnames(df_pbmc))) df_pbmc$gene <- rownames(df_pbmc)
if (!("gene" %in% colnames(df_vcc)))  df_vcc$gene  <- rownames(df_vcc)

df <- merge(
  df_pbmc[, c("gene", "avg_log2FC", "p_val_adj")],
  df_vcc[,  c("gene", "avg_log2FC", "p_val_adj")],
  by = "gene", suffixes = c("_PBMC", "_VCC")
)

df$p_val_adj_PBMC[is.na(df$p_val_adj_PBMC)] <- 1
df$p_val_adj_VCC[is.na(df$p_val_adj_VCC)]   <- 1

# ==============================================================================
# Classify
# ==============================================================================

df$sig_PBMC <- df$p_val_adj_PBMC < PVAL_CUTOFF
df$sig_VCC  <- df$p_val_adj_VCC  < PVAL_CUTOFF

df$category <- with(df,
  ifelse(sig_PBMC & sig_VCC,   "Significant in both",
  ifelse(sig_PBMC & !sig_VCC,  "Significant in PBMC only",
  ifelse(!sig_PBMC & sig_VCC,  "Significant in VCC only",
         "Not significant in either"))))

df$category <- factor(df$category, levels = c(
  "Significant in both",
  "Significant in PBMC only",
  "Significant in VCC only",
  "Not significant in either"
))

# Label all significant genes (ggrepel suppresses overlapping ones)
df$label <- ifelse(df$category != "Not significant in either", df$gene, NA_character_)

message(sprintf("Genes after merge: %d", nrow(df)))
print(table(df$category))

# ==============================================================================
# Plot
# ==============================================================================

p <- ggplot(df, aes(x = avg_log2FC_PBMC, y = avg_log2FC_VCC, color = category)) +
  geom_point(alpha = 0.75, size = 1.6, shape = 16) +
  geom_vline(xintercept = c(-FC_LINE, FC_LINE),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-FC_LINE, FC_LINE),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c(
    "Significant in both"         = "red3",
    "Significant in PBMC only"    = "deepskyblue4",
    "Significant in VCC only"     = "orange3",
    "Not significant in either"   = "grey"
  )) +
  geom_text_repel(
    aes(label = label),
    size = 6, force = 0.1, max.overlaps = 8,
    box.padding = 0.2, point.padding = 0.15, segment.size = 0.2,
    na.rm = TRUE, color = "black", segment.color = "black"
  ) +
  labs(x = "Average log2FC (PBMCs)", y = "Average log2FC (VCC)") +
  theme_classic() +
  theme(
    plot.title        = element_blank(),
    axis.title        = element_text(size = 20, color = "black"),
    axis.text         = element_text(size = 16, color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "none"
  )

ggsave(file.path(FIGURES_DIR, "fig1i_scatter_BC002vsBC011_VCC_vs_PBMC.pdf"),
       p, width = 5, height = 5, device = cairo_pdf)
message("Saved: fig1i_scatter_BC002vsBC011_VCC_vs_PBMC.pdf")

message("Done.")
print(sessionInfo())
