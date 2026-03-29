#!/usr/bin/env Rscript
# ==============================================================================
# fig1d_violin_top_genes.R
# Violin plots of expression by probe barcode for SRP9 and SPCS1 (Figure 1D)
# ==============================================================================
#
# Loads the eta-squared results to confirm rank order, then plots per-cell
# expression of the top two genes (SRP9, rank 1; SPCS1, rank 2) across probe
# barcodes in NTC cells. Violins are split and colored by sequencing lane.
#
# Uses Seurat's VlnPlot() on the NTC-subset Seurat object.
#
# INPUT
#   results/seurat_vcc_annotated.qs              — from 01_setup_vcc_seurat.R
#   results/variance_explained_probe_barcode_NTC.csv  — from 04_compute_variance_explained.R
#
# OUTPUT
#   results/figures/fig1d_vln_SRP9_by_probe_barcode.pdf
#   results/figures/fig1d_vln_SPCS1_by_probe_barcode.pdf
#
# USAGE
#   Rscript scripts/figure1/plot/fig1d_violin_top_genes.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(ggthemes)
  library(dplyr)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

VCC_QS  <- file.path(RESULTS_DIR, "seurat_vcc_annotated.qs")
ETA_CSV <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_NTC.csv")
NTC_LABEL <- "non-targeting"

# ==============================================================================
# Confirm top genes from eta-squared results
# ==============================================================================

if (!file.exists(ETA_CSV)) {
  stop("Eta-squared CSV not found: ", ETA_CSV,
       "\nRun 04_compute_variance_explained.R first.")
}

df_eta    <- read.csv(ETA_CSV) %>% arrange(rank_bc)
top2      <- df_eta$gene[1:2]
message(sprintf("Top 2 genes by probe barcode eta²: %s", paste(top2, collapse = ", ")))

# ==============================================================================
# Load Seurat and subset to NTC cells
# ==============================================================================

if (!file.exists(VCC_QS)) {
  stop("VCC Seurat not found: ", VCC_QS, "\nRun 01_setup_vcc_seurat.R first.")
}

message("Loading VCC Seurat object...")
seurat_vcc <- qread(VCC_QS)

cells_ntc  <- colnames(seurat_vcc)[
  !is.na(seurat_vcc$target_gene) &
  seurat_vcc$target_gene == NTC_LABEL &
  !is.na(seurat_vcc$probe_barcode)
]
seurat_ntc <- subset(seurat_vcc, cells = cells_ntc)
rm(seurat_vcc); gc()

message(sprintf("NTC cells: %d", ncol(seurat_ntc)))

Idents(seurat_ntc) <- seurat_ntc$probe_barcode

# Lane colors (colorblind-friendly palette)
lane_colors <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73")
seurat_ntc$lane <- as.factor(seurat_ntc$lane)

vln_theme <- theme_few() +
  theme(
    plot.title        = element_text(size = 15, color = "black", hjust = 0.5,
                                     face = "italic"),
    axis.title        = element_text(size = 17, color = "black"),
    axis.text         = element_text(size = 14, color = "black"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.title      = element_text(size = 13, color = "black"),
    legend.text       = element_text(size = 13, color = "black")
  )

# ==============================================================================
# Plot and save
# ==============================================================================

for (gene in top2) {
  if (!(gene %in% rownames(seurat_ntc))) {
    message(sprintf("Skipping %s: not found in matrix", gene))
    next
  }
  eta_val  <- df_eta$eta_sq_bc[df_eta$gene == gene]
  rank_val <- df_eta$rank_bc[df_eta$gene == gene]

  p <- VlnPlot(seurat_ntc, features = gene, split.by = "lane", pt.size = 0) +
    labs(
      title = sprintf("%s  (rank %d, \u03b7\u00b2 = %.3f)", gene, rank_val, eta_val),
      x     = "Probe barcode",
      y     = "Normalized expression (NTC cells)",
      fill  = "Lane"
    ) +
    scale_fill_manual(values  = lane_colors) +
    scale_color_manual(values = lane_colors) +
    guides(color = guide_legend(title = "Lane"),
           fill  = guide_legend(title = "Lane")) +
    vln_theme

  out_pdf <- file.path(FIGURES_DIR,
    sprintf("fig1d_vln_%s_by_probe_barcode.pdf", gene))
  ggsave(out_pdf, plot = p, width = 12, height = 3.5, device = cairo_pdf)
  message("Saved: ", basename(out_pdf))
}

message("Done.")
print(sessionInfo())
