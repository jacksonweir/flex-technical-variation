#!/usr/bin/env Rscript
# ==============================================================================
# fig1bc_volcanos.R
# Volcano plots — lane effect vs probe barcode effect (Figures 1B and 1C)
# ==============================================================================
#
# Generates two volcano plots with shared axis limits for direct comparison:
#
#   Figure 1B — Lane 1 vs Lane 2, BC001 cells only
#     Represents the lane-to-lane technical variation within a single barcode.
#
#   Figure 1C — BC002 vs BC011, Lane 1 cells only
#     Represents the probe barcode technical variation within a single lane.
#
# Significance threshold: p_val_adj < 0.05 and |log2FC| >= 0.3.
# Points colored: tomato3 (up), steelblue3 (down), grey (not significant).
# Significant genes labelled with ggrepel.
#
# INPUT
#   results/DE_lane1_vs_lane2_BC001_wilcox.csv  — from 02_compute_pairwise_de.R
#   results/DE_BC002_vs_BC011_lane1_wilcox.csv  — from 02_compute_pairwise_de.R
#
# OUTPUT
#   results/figures/fig1b_volcano_lane1vsLane2_BC001only.pdf
#   results/figures/fig1c_volcano_BC002vsBC011_lane1only.pdf
#
# USAGE
#   Rscript scripts/figure1/plot/fig1bc_volcanos.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
})

RESULTS_DIR  <- "results"
FIGURES_DIR  <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

PVAL_CUTOFF <- 0.05
FC_THRESH   <- 0.3

# ==============================================================================
# Load DE results
# ==============================================================================

de_lane <- read.csv(file.path(RESULTS_DIR, "DE_lane1_vs_lane2_BC001_wilcox.csv"))
de_bc   <- read.csv(file.path(RESULTS_DIR, "DE_BC002_vs_BC011_lane1_wilcox.csv"))

message(sprintf("Lane DE: %d genes  |  BC DE: %d genes", nrow(de_lane), nrow(de_bc)))

# ==============================================================================
# Annotate significance direction
# ==============================================================================

annotate_de <- function(de, pval_cutoff, fc_thresh) {
  de$log10_pval <- -log10(de$p_val_adj + 1e-300)
  de$reg_dir <- "Not Significant"
  de$reg_dir[de$p_val_adj < pval_cutoff & de$avg_log2FC >  fc_thresh] <- "Upregulated"
  de$reg_dir[de$p_val_adj < pval_cutoff & de$avg_log2FC < -fc_thresh] <- "Downregulated"
  de$label <- ifelse(de$reg_dir != "Not Significant", de$gene, NA_character_)
  de
}

de_lane <- annotate_de(de_lane, PVAL_CUTOFF, FC_THRESH)
de_bc   <- annotate_de(de_bc,   PVAL_CUTOFF, FC_THRESH)

for (nm in c("Lane", "BC")) {
  de <- if (nm == "Lane") de_lane else de_bc
  message(sprintf("%s — Up: %d  Down: %d  NS: %d",
    nm, sum(de$reg_dir == "Upregulated"),
    sum(de$reg_dir == "Downregulated"),
    sum(de$reg_dir == "Not Significant")))
}

# ==============================================================================
# Shared axis limits (driven by the more extreme dataset)
# ==============================================================================

all_lfc  <- c(de_lane$avg_log2FC, de_bc$avg_log2FC)
all_yval <- c(de_lane$log10_pval, de_bc$log10_pval)

xlims <- c(floor(min(all_lfc)   * 10) / 10,
           ceiling(max(all_lfc) * 10) / 10)
ylims <- c(0, ceiling(max(all_yval) / 10) * 10)

message(sprintf("Shared xlim: [%.1f, %.1f]  ylim: [0, %d]", xlims[1], xlims[2], ylims[2]))

# ==============================================================================
# Volcano plot function
# ==============================================================================

make_volcano <- function(de, x_label, xlims, ylims) {
  de <- de[order(de$reg_dir == "Not Significant", decreasing = TRUE), ]

  ggplot(de, aes(x = avg_log2FC, y = log10_pval)) +
    geom_point(aes(color = reg_dir), alpha = 0.8, size = 1, shape = 16) +
    scale_color_manual(values = c(
      "Upregulated"     = "tomato3",
      "Downregulated"   = "steelblue3",
      "Not Significant" = "grey"
    )) +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(PVAL_CUTOFF),
               linetype = "dashed", color = "black") +
    geom_text_repel(
      aes(label = label),
      size = 3, force = 0.5, max.overlaps = 8,
      box.padding = 0.1, point.padding = 0.1, segment.size = 0.1,
      na.rm = TRUE, color = "black", segment.color = "black"
    ) +
    scale_x_continuous(limits = xlims) +
    scale_y_continuous(limits = ylims) +
    labs(x = x_label, y = "-log10(adjusted p-value)") +
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
}

# ==============================================================================
# Save
# ==============================================================================

p_b <- make_volcano(de_lane,
  x_label = "Average log2 Fold Change (Lane 1 / Lane 2)",
  xlims = xlims, ylims = ylims)
ggsave(file.path(FIGURES_DIR, "fig1b_volcano_lane1vsLane2_BC001only.pdf"),
       p_b, width = 7, height = 5, device = cairo_pdf)
message("Saved: fig1b_volcano_lane1vsLane2_BC001only.pdf")

p_c <- make_volcano(de_bc,
  x_label = "Average log2 Fold Change (BC002 / BC011)",
  xlims = xlims, ylims = ylims)
ggsave(file.path(FIGURES_DIR, "fig1c_volcano_BC002vsBC011_lane1only.pdf"),
       p_c, width = 7, height = 5, device = cairo_pdf)
message("Saved: fig1c_volcano_BC002vsBC011_lane1only.pdf")

message("Done.")
print(sessionInfo())
