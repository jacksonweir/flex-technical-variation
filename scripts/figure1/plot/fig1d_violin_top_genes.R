#!/usr/bin/env Rscript
# ==============================================================================
# fig1d_violin_top_genes.R
# Violin plots of expression by probe barcode for SRP9 and SPCS1 (Figure 1D)
# ==============================================================================
#
# Plots per-cell expression of the top two η²-ranked genes (SRP9, rank 1;
# SPCS1, rank 2) across probe barcodes in NTC cells. Violins are split and
# colored by sequencing lane.
#
# Reads from the pre-extracted CSV (no Seurat required at plot time). Run
# scripts/figure1/00_prep_violin_data.R first to generate the input CSV.
#
# INPUT
#   results/vcc_ntc_violin_data.csv              — from 00_prep_violin_data.R
#   results/variance_explained_probe_barcode_NTC.csv  — from 04_compute_variance_explained.R
#
# OUTPUT
#   results/figures/fig1d_SRP9.pdf
#   results/figures/fig1d_SPCS1.pdf
#
# USAGE
#   Rscript scripts/figure1/plot/fig1d_violin_top_genes.R
#
# PACKAGES
#   ggplot2, dplyr
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

VIOLIN_CSV <- file.path(RESULTS_DIR, "vcc_ntc_violin_data.csv")
ETA_CSV    <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_NTC.csv")

for (f in c(VIOLIN_CSV, ETA_CSV)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f,
         "\nRun 00_prep_violin_data.R and 04_compute_variance_explained.R first.")
  }
}

# ==============================================================================
# Load data
# ==============================================================================

df  <- read.csv(VIOLIN_CSV, stringsAsFactors = FALSE)
eta <- read.csv(ETA_CSV) |> dplyr::arrange(rank_bc)

top_genes <- unique(df$gene)
message("Genes to plot: ", paste(top_genes, collapse = ", "))

# ==============================================================================
# Factor levels and colors
# ==============================================================================

bc_order         <- paste0("BC", sprintf("%03d", 1:16))
df$probe_barcode <- factor(df$probe_barcode, levels = bc_order)
df$lane          <- factor(df$lane, levels = c("1", "2", "3"))

lane_colors <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73")

# ==============================================================================
# Theme
# ==============================================================================

vln_theme <- theme_classic() +
  theme(
    plot.title        = element_text(size = 30, color = "black", hjust = 0.5,
                                     face = "italic"),
    axis.title        = element_text(size = 28, color = "black"),
    axis.text         = element_text(size = 24, color = "black"),
    axis.text.x       = element_text(size = 24, angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.text       = element_text(size = 24, color = "black"),
    legend.title      = element_text(size = 26, color = "black")
  )

# ==============================================================================
# Plot and save (one PDF per gene)
# ==============================================================================

for (gene in top_genes) {
  eta_val  <- eta$eta_sq_bc[eta$gene == gene]
  rank_val <- eta$rank_bc[eta$gene == gene]
  if (length(eta_val) == 0) {
    message("Gene not found in \u03b7\u00b2 table: ", gene, " \u2014 skipping")
    next
  }

  df_gene <- df[df$gene == gene & !is.na(df$probe_barcode), ]

  # y-axis breaks every 1 unit
  y_range  <- range(df_gene$expression, na.rm = TRUE)
  y_breaks <- seq(floor(y_range[1]), ceiling(y_range[2]), by = 1)
  y_labels <- sprintf("%.1f", y_breaks)

  p <- ggplot(df_gene, aes(x = probe_barcode, y = expression, fill = lane)) +
    geom_violin(color = "black", scale = "width", trim = TRUE,
                linewidth = 0.3, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = lane_colors, name = "Lane") +
    scale_y_continuous(breaks = y_breaks, labels = y_labels) +
    labs(
      title = gene,
      x     = "Probe barcode",
      y     = "Expression"
    ) +
    vln_theme

  out_pdf <- file.path(FIGURES_DIR, sprintf("fig1d_%s.pdf", gene))
  ggsave(out_pdf, plot = p, width = 12, height = 3.5, device = cairo_pdf)
  message("Saved: ", basename(out_pdf))
}

message("Done.")
print(sessionInfo())
