#!/usr/bin/env Rscript
# ==============================================================================
# fig2c_violin_top_genes.R
# Violin plots of SRP9 and SPCS1 by probe barcode — A375 FLEXv2 (Figure 2C)
# ==============================================================================
#
# Plots expression distributions of SRP9 and SPCS1 across Flex v2 probe
# barcodes in A375 cells. These genes show substantial barcode-driven expression
# variation in the VCC Flex v1 dataset (Figure 1D). Plotting them in the A375
# Flex v2 dataset allows direct comparison of the technical artifact magnitude
# between kit versions.
#
# Reads from the pre-extracted CSV (no Seurat required at plot time). Run
# scripts/figure2/00_prep_violin_data_flexv2.R first (mel_spatial env) to
# generate the input CSV.
#
# INPUT
#   results/a375_flexv2_violin_data.csv   — from 00_prep_violin_data_flexv2.R
#
# OUTPUT
#   results/figures/fig2c_SRP9.pdf
#   results/figures/fig2c_SPCS1.pdf
#
# USAGE
#   conda run -n trekker Rscript scripts/figure2/plot/fig2c_violin_top_genes.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

VIOLIN_CSV <- file.path(RESULTS_DIR, "a375_flexv2_violin_data.csv")

if (!file.exists(VIOLIN_CSV)) {
  stop("Violin data not found: ", VIOLIN_CSV,
       "\nRun scripts/figure2/00_prep_violin_data_flexv2.R first (mel_spatial env).")
}

# ==============================================================================
# Load data
# ==============================================================================

df <- read.csv(VIOLIN_CSV, stringsAsFactors = FALSE)

top_genes <- unique(df$gene)
message("Genes to plot: ", paste(top_genes, collapse = ", "))

# ==============================================================================
# Factor levels
# ==============================================================================

bc_order         <- c(paste0("A-A", sprintf("%02d", 1:12)),
                      paste0("A-B", sprintf("%02d", 1:4)))
df$probe_barcode <- factor(df$probe_barcode, levels = bc_order)

# ==============================================================================
# Theme
# ==============================================================================

vln_theme <- theme_classic() +
  theme(
    plot.title        = element_text(size = 34, color = "black", hjust = 0.5,
                                     face = "italic"),
    axis.title        = element_text(size = 30, color = "black"),
    axis.text         = element_text(size = 26, color = "black"),
    axis.text.x       = element_text(size = 26, angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "none"
  )

# ==============================================================================
# Plot and save (one PDF per gene)
# ==============================================================================

for (gene in top_genes) {
  df_gene <- df[df$gene == gene & !is.na(df$probe_barcode), ]

  # y-axis breaks every 1 unit, 1 decimal place labels
  y_range  <- range(df_gene$expression, na.rm = TRUE)
  y_breaks <- seq(floor(y_range[1]), ceiling(y_range[2]), by = 1)
  y_labels <- sprintf("%.1f", y_breaks)

  p <- ggplot(df_gene, aes(x = probe_barcode, y = expression)) +
    geom_violin(fill = "steelblue3", color = "black", scale = "width",
                trim = TRUE, linewidth = 0.3) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels,
                       minor_breaks = NULL) +
    labs(
      title = gene,
      x     = "Probe barcode",
      y     = "Expression"
    ) +
    vln_theme

  out_pdf <- file.path(FIGURES_DIR, sprintf("fig2c_%s.pdf", gene))
  ggsave(out_pdf, plot = p, width = 12, height = 3.5, device = cairo_pdf)
  message("Saved: ", basename(out_pdf))
}

message("Done.")
print(sessionInfo())
