#!/usr/bin/env Rscript
# ==============================================================================
# fig2c_violin_top_genes.R
# Violin plots of SRP9 and SPCS1 by probe barcode — A375 FLEXv2 (Figure 2C)
# ==============================================================================
#
# Plots expression distributions of SRP9 and SPCS1 across Flex v2 probe barcode
# numbers in PBS-treated A375 cells. These genes show substantial barcode-driven
# expression variation in the VCC Flex v1 dataset (Figure 1D). Plotting them in
# a Flex v2 PBS-cell dataset allows direct comparison of the technical artifact
# magnitude between probe kit versions.
#
# Cells are from a controlled experiment in which A375 cells were treated with
# PBS and distributed across multiple Flex v2 sample barcodes. Because all cells
# received the same treatment, any expression variation across probe barcodes is
# a technical artifact of the barcoding chemistry.
#
# NOTE: This script requires the PBS-cell Seurat object from the Chen lab
# servers. Set A375_PBS_RDS below to the correct path. This object uses
# `probe_barcode_number` as the grouping variable (Flex v2 barcode labels).
#
# INPUT  (user must set path)
#   A375_PBS_RDS — RDS path to flex_v2_seurat_obj_filt_PBS.rds
#
# OUTPUT
#   results/figures/fig2c_vln_SRP9_by_probe_barcode_a375.pdf
#   results/figures/fig2c_vln_SPCS1_by_probe_barcode_a375.pdf
#
# USAGE
#   Rscript scripts/figure2/plot/fig2c_violin_top_genes.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggthemes)
})

# Set to path of the PBS-cell A375 FLEXv2 Seurat object (Chen lab servers)
A375_PBS_RDS <- "/path/to/flex_v2_seurat_obj_filt_PBS.rds"

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

GENES <- c("SRP9", "SPCS1")

if (!file.exists(A375_PBS_RDS)) {
  stop("A375 PBS Seurat object not found: ", A375_PBS_RDS,
       "\nThis data is stored on the Chen lab servers (not publicly available).",
       "\nUpdate A375_PBS_RDS at the top of this script.")
}

# ==============================================================================
# Load and plot
# ==============================================================================

message("Loading A375 FLEXv2 PBS-cell Seurat object...")
obj <- readRDS(A375_PBS_RDS)
message(sprintf("Cells: %d", ncol(obj)))

genes_present <- GENES[GENES %in% rownames(obj)]
if (length(genes_present) == 0) stop("None of the requested genes found in object.")
if (length(genes_present) < length(GENES)) {
  message("Note: genes not found in object: ",
          paste(setdiff(GENES, genes_present), collapse = ", "))
}

Idents(obj) <- obj$probe_barcode_number
bc_vec  <- factor(obj$probe_barcode_number,
                  levels = sort(unique(obj$probe_barcode_number)))
n_bc    <- nlevels(bc_vec)
bc_cols <- rep("steelblue3", n_bc)

vln_theme <- theme_few() +
  theme(
    plot.title        = element_text(size = 17, color = "black", hjust = 0.5,
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
    legend.position   = "none"
  )

for (gene in genes_present) {
  message(sprintf("Plotting: %s", gene))
  p <- VlnPlot(obj, features = gene, pt.size = 0, cols = bc_cols) +
    labs(title = gene, x = "Probe barcode", y = "Expression") +
    vln_theme

  out_file <- file.path(FIGURES_DIR,
    sprintf("fig2c_vln_%s_by_probe_barcode_a375.pdf", gene))
  ggsave(out_file, plot = p, width = 10, height = 3, device = cairo_pdf)
  message("Saved: ", basename(out_file))
}

rm(obj); gc()
message("Done.")
print(sessionInfo())
