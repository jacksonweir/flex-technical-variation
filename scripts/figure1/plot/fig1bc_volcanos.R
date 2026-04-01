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
# Significant genes labelled with ggrepel (top 60 by score; SRP9/ICE2/SPCS1 forced).
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
# PACKAGES
#   ggplot2, ggthemes, ggrepel
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
MAX_LABELS  <- 60

# ==============================================================================
# Load DE results
# ==============================================================================

de_lane <- read.csv(file.path(RESULTS_DIR, "DE_lane1_vs_lane2_BC001_wilcox.csv"))
de_bc   <- read.csv(file.path(RESULTS_DIR, "DE_BC002_vs_BC011_lane1_wilcox.csv"))

message(sprintf("Lane DE: %d genes  |  BC DE: %d genes", nrow(de_lane), nrow(de_bc)))

# ==============================================================================
# Annotate significance direction and select labels
# ==============================================================================

annotate_de <- function(de) {
  de$log10_pval <- -log10(de$p_val_adj + 1e-300)
  de$reg_dir    <- "Not Significant"
  de$reg_dir[de$p_val_adj < PVAL_CUTOFF & de$avg_log2FC >  FC_THRESH] <- "Upregulated"
  de$reg_dir[de$p_val_adj < PVAL_CUTOFF & de$avg_log2FC < -FC_THRESH] <- "Downregulated"

  # Label top genes by combined score (p-value × effect size)
  sig <- de$reg_dir != "Not Significant"
  if (sum(sig) > MAX_LABELS) {
    score     <- de$log10_pval * abs(de$avg_log2FC)
    top_genes <- de$gene[sig][order(score[sig], decreasing = TRUE)[seq_len(MAX_LABELS)]]
    de$label  <- ifelse(de$gene %in% top_genes, de$gene, NA_character_)
  } else {
    de$label <- ifelse(sig, de$gene, NA_character_)
  }
  de
}

de_lane <- annotate_de(de_lane)
de_bc   <- annotate_de(de_bc)
# Suppress two crowded labels in Fig 1C (will be nearby labeled genes)
de_bc$label[de_bc$gene %in% c("COMMD7", "CSNK1E")] <- NA_character_

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

xlims <- c(floor(min(all_lfc,  na.rm = TRUE) * 10) / 10,
           ceiling(max(all_lfc, na.rm = TRUE) * 10) / 10)
ylims <- c(0, ceiling(max(all_yval, na.rm = TRUE) / 10) * 10)

message(sprintf("Shared xlim: [%.1f, %.1f]  ylim: [0, %d]", xlims[1], xlims[2], ylims[2]))

# ==============================================================================
# Volcano plot function
# ==============================================================================

make_volcano <- function(de, x_label, xlims, ylims, force_labels = character(0)) {
  de <- de[order(de$reg_dir == "Not Significant", decreasing = TRUE), ]

  # Separate forced labels to guarantee they appear regardless of overlap limit
  de$label_regular <- ifelse(de$gene %in% force_labels, NA_character_, de$label)
  de_forced        <- de[de$gene %in% force_labels, ]

  ggplot(de, aes(x = avg_log2FC, y = log10_pval)) +
    geom_point(aes(color = reg_dir), alpha = 0.8, size = 2, shape = 16) +
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
      aes(label = label_regular),
      size = 6, force = 2.5, max.overlaps = 5,
      box.padding = 0.3, point.padding = 0.2, segment.size = 0.2,
      min.segment.length = 0, na.rm = TRUE,
      color = "black", segment.color = "black"
    ) +
    geom_text_repel(
      data = de_forced,
      aes(label = gene),
      size = 6, force = 2.5, max.overlaps = Inf,
      box.padding = 0.3, point.padding = 0.2, segment.size = 0.2,
      min.segment.length = 0,
      color = "black", segment.color = "black",
      seed = 42
    ) +
    scale_x_continuous(limits = xlims) +
    scale_y_continuous(limits = ylims) +
    labs(x = x_label, y = "-log10(adjusted p-value)") +
    theme_few() +
    theme(
      axis.title        = element_text(size = 20, color = "black"),
      axis.text         = element_text(size = 16, color = "black"),
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
  xlims = xlims, ylims = ylims,
  force_labels = c("SRP9", "ICE2", "SPCS1"))
ggsave(file.path(FIGURES_DIR, "fig1c_volcano_BC002vsBC011_lane1only.pdf"),
       p_c, width = 7, height = 5, device = cairo_pdf)
message("Saved: fig1c_volcano_BC002vsBC011_lane1only.pdf")

message("Done.")
print(sessionInfo())
