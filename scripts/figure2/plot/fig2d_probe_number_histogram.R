#!/usr/bin/env Rscript
# ==============================================================================
# fig2d_probe_number_histogram.R
# Probe count per gene — Flex v1 vs Flex v2 side-by-side histograms (Figure 2D)
# ==============================================================================
#
# For each expressed gene, counts the number of unique probes targeting it in
# the Flex v1 (VCC) and Flex v2 (A375) probe sets, then plots side-by-side
# histograms stratified by DE status. This tests whether genes targeted by more
# probes are less likely to show probe-barcode-associated technical variation.
#
# DE gene classification:
#   Flex v1 (VCC):  p_val_adj < 0.05 AND |avg_log2FC| >= 0.3 from a
#     FindAllMarkers across all 16 probe barcodes in NTC cells (all lanes).
#     Input: VCC_WILCOX_CSV — see note below.
#   Flex v2 (A375): p_val_adj < 0.05 AND |avg_log2FC| >= 0.3 from a
#     FindAllMarkers across Flex v2 sample barcodes in PBS-treated A375 cells.
#     Input: A375_WILCOX_CSV — see note below.
#
# Expressed genes filter:
#   Flex v1 (VCC):  genes expressed in >= 10% of NTC cells
#   Flex v2 (A375): genes with pct_expressed >= 0.10 in the PBS-cell dataset
#                   (derived from results/variance_explained_probe_barcode_A375_flexv2.csv)
#
# NOTE ON DATA ACCESS: VCC_WILCOX_CSV and A375_WILCOX_CSV are pre-computed DE
# result files generated from Chen lab analyses on server data. They are not
# reproduced by other scripts in this repository. Update the path constants
# below if running from a different location. The VCC_MIN01_CSV file is derived
# from the VCC NTC cells (genes expressed >= 10%); update accordingly.
#
# Probe manifests must be downloaded (see data/README.md):
#   Flex v1: 10x Chromium Human Transcriptome Probe Set v1.0.1
#   Flex v2: 10x Chromium Human Transcriptome Probe Set v2.0.0
#
# INPUT
#   data/vcc_probe_manifest.csv          — Flex v1 probe manifest (public, 10x)
#   data/a375_v2_probe_manifest.csv      — Flex v2 probe manifest (public, 10x)
#   VCC_MIN01_CSV   — genes expressed >= 10% of NTC cells (Chen lab)
#   VCC_WILCOX_CSV  — VCC FindAllMarkers by probe barcode, NTC cells (Chen lab)
#   A375_WILCOX_CSV — A375 FindAllMarkers by Flex v2 barcode, PBS cells (Chen lab)
#   results/variance_explained_probe_barcode_A375_flexv2.csv  — from 10_compute...
#
# OUTPUT
#   results/figures/fig2d_probe_number_histogram_v1_vs_v2.pdf
#
# USAGE
#   Rscript scripts/figure2/plot/fig2d_probe_number_histogram.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Paths to Chen lab reference files — update if running outside Chen lab servers
VCC_MIN01_CSV   <- "/path/to/vcc_genes_min_0.1.csv"
VCC_WILCOX_CSV  <- "/path/to/260302_VCC_wilcox_de_genes.csv"
A375_WILCOX_CSV <- "/path/to/260222_flex_v2_PBS_de_wilcox.csv"

# Probe manifests — download per data/README.md and place at these paths
VCC_MANIFEST_CSV  <- "data/vcc_probe_manifest.csv"
A375_MANIFEST_CSV <- "data/a375_v2_probe_manifest.csv"

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

A375_ETA_CSV  <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_A375_flexv2.csv")
LOG2FC_CUTOFF <- 0.3
PVAL_CUTOFF   <- 0.05

for (f in c(VCC_MIN01_CSV, VCC_WILCOX_CSV, A375_WILCOX_CSV,
            VCC_MANIFEST_CSV, A375_MANIFEST_CSV, A375_ETA_CSV)) {
  if (!file.exists(f)) stop("Required file not found: ", f)
}

# ==============================================================================
# VCC (Flex v1)
# ==============================================================================

message("Processing VCC (Flex v1) data...")

manifest_vcc <- read.csv(VCC_MANIFEST_CSV, comment.char = "#")
manifest_vcc <- manifest_vcc[manifest_vcc$included == TRUE, ]
manifest_vcc$gene_name <- sapply(strsplit(manifest_vcc$probe_id, "\\|"), `[`, 2)

probe_counts_vcc <- manifest_vcc %>%
  group_by(gene_name) %>%
  summarise(n_probe_ids = n_distinct(probe_id), .groups = "drop")

min01_genes_vcc  <- read.csv(VCC_MIN01_CSV)$gene
probe_counts_vcc <- probe_counts_vcc[
  probe_counts_vcc$gene_name %in% min01_genes_vcc, ]

wilcox_vcc      <- read.csv(VCC_WILCOX_CSV)
de_genes_vcc    <- unique(wilcox_vcc$gene[
  !is.na(wilcox_vcc$p_val_adj) & wilcox_vcc$p_val_adj < PVAL_CUTOFF &
  abs(wilcox_vcc$avg_log2FC) > LOG2FC_CUTOFF])

probe_counts_vcc$is_de <- ifelse(
  probe_counts_vcc$gene_name %in% de_genes_vcc, "DE", "Not DE")

message(sprintf("VCC: %d expressed genes | DE: %d | Not DE: %d",
  nrow(probe_counts_vcc),
  sum(probe_counts_vcc$is_de == "DE"),
  sum(probe_counts_vcc$is_de == "Not DE")))

# ==============================================================================
# A375 (Flex v2)
# ==============================================================================

message("Processing A375 (Flex v2) data...")

manifest_a375 <- read.csv(A375_MANIFEST_CSV, comment.char = "#")
manifest_a375 <- manifest_a375[manifest_a375$included == TRUE, ]

probe_counts_a375 <- manifest_a375 %>%
  group_by(gene_name) %>%
  summarise(n_probe_ids = n_distinct(probe_id), .groups = "drop")

a375_eta         <- read.csv(A375_ETA_CSV)
min01_genes_a375 <- a375_eta$gene[a375_eta$pct_expressed >= 0.10]
probe_counts_a375 <- probe_counts_a375[
  probe_counts_a375$gene_name %in% min01_genes_a375, ]

wilcox_a375   <- read.csv(A375_WILCOX_CSV)
de_genes_a375 <- unique(wilcox_a375$gene[
  !is.na(wilcox_a375$p_val_adj) & wilcox_a375$p_val_adj < PVAL_CUTOFF &
  abs(wilcox_a375$avg_log2FC) > LOG2FC_CUTOFF])

probe_counts_a375$is_de <- ifelse(
  probe_counts_a375$gene_name %in% de_genes_a375, "DE", "Not DE")

message(sprintf("A375: %d expressed genes | DE: %d | Not DE: %d",
  nrow(probe_counts_a375),
  sum(probe_counts_a375$is_de == "DE"),
  sum(probe_counts_a375$is_de == "Not DE")))

# ==============================================================================
# Plot
# ==============================================================================

make_panel <- function(probe_counts, panel_title) {
  probe_counts$is_de <- factor(probe_counts$is_de, levels = c("Not DE", "DE"))
  ggplot(probe_counts, aes(x = n_probe_ids, fill = is_de)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                   color = "black", linewidth = 0.3) +
    facet_wrap(~ is_de, ncol = 1) +
    scale_fill_manual(values = c("Not DE" = "grey", "DE" = "tomato3"),
                      name = NULL) +
    labs(x = "Number of probes per gene", y = "Density",
         title = panel_title) +
    theme_classic() +
    theme(
      plot.title           = element_text(size = 18, hjust = 0.5),
      axis.title           = element_text(size = 16),
      axis.text            = element_text(size = 14),
      strip.text           = element_blank(),
      strip.background     = element_blank(),
      legend.position      = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.key.spacing.y = unit(4, "pt"),
      legend.text          = element_text(size = 14)
    )
}

p <- make_panel(probe_counts_vcc,  "Flex v1") +
     make_panel(probe_counts_a375, "Flex v2") +
     plot_layout(guides = "collect") &
     theme(legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "fig2d_probe_number_histogram_v1_vs_v2.pdf"),
       p, width = 10, height = 5, device = cairo_pdf)
message("Saved: fig2d_probe_number_histogram_v1_vs_v2.pdf")

message("Done.")
print(sessionInfo())
