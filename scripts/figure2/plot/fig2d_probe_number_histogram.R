#!/usr/bin/env Rscript
# ==============================================================================
# fig2d_probe_number_histogram.R
# Probe count per gene — Flex v1 vs Flex v2 side-by-side histograms (Figure 2D)
# ==============================================================================
#
# For each expressed gene, counts the number of unique probes targeting it in
# the Flex v1 (VCC) and Flex v2 (A375) probe sets, then plots histograms
# stratified by DE status (DE = differentially expressed across probe barcodes;
# Not DE = not). This tests whether genes with more probes are less likely to
# show probe-barcode-associated technical variation.
#
# DE gene classification per dataset:
#   Flex v1 (VCC): p_val_adj < 0.05 AND |avg_log2FC| >= 0.3 in at least one
#     probe barcode comparison from FindAllMarkers (Lane 1, NTC cells only)
#   Flex v2 (A375): p_val_adj < 0.05 AND |avg_log2FC| >= 0.3 in at least one
#     probe barcode comparison from FindAllMarkers
#
# NOTE ON DISCREPANCY: In the published figure, the A375 (Flex v2) DE gene
# classification was performed on PBS-only cells (same treatment, multiple Flex
# v2 probe barcodes) to match the controlled VCC NTC design. The script here
# uses all A375 cells grouped by ligand condition (from 09_compute_a375_de.R),
# which conflates biological and technical effects. For an analysis that exactly
# matches the paper, the A375 DE gene list should be derived from PBS-only cells
# using the colleague's analysis (data available on Chen lab servers only).
#
# Probe manifests must be downloaded by the user (see data/README.md):
#   Flex v1 probe set: 10x Genomics Chromium Human Transcriptome Probe Set v1.0.1
#   Flex v2 probe set: 10x Genomics Chromium Human Transcriptome Probe Set v2.0.0
#
# INPUT
#   results/variance_explained_probe_barcode_NTC.csv    — from 04_compute_variance_explained.R
#   results/DE_FindAllMarkers_lane1_wilcox.csv           — from 03_compute_findallmarkers_lane1.R
#   results/variance_explained_probe_barcode_A375_flexv2.csv  — from 10_compute_variance_explained_a375.R
#   results/DE_A375_flexv2_FindAllMarkers_wilcox.csv     — from 09_compute_a375_de.R
#   VCC_MANIFEST_CSV  — Flex v1 probe manifest (user must download, see below)
#   A375_MANIFEST_CSV — Flex v2 probe manifest (user must download, see below)
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

# Paths to probe manifests — update these after downloading (see data/README.md)
# Flex v1: Chromium Human Transcriptome Probe Set v1.0.1 (GRCh38-2020-A)
#   Download from 10x Genomics support site
VCC_MANIFEST_CSV  <- "data/vcc_probe_manifest.csv"
# Flex v2: Chromium Human Transcriptome Probe Set v2.0.0 (GRCh38-2024-A)
#   Download from: https://www.10xgenomics.com/support/software/cell-ranger/downloads
A375_MANIFEST_CSV <- "data/a375_v2_probe_manifest.csv"

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

LOG2FC_CUTOFF <- 0.3
PVAL_CUTOFF   <- 0.05
MIN_PCT_EXPR  <- 0.10   # minimum fraction of cells expressing the gene

# Input CSVs
VCC_ETA_CSV  <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_NTC.csv")
VCC_DE_CSV   <- file.path(RESULTS_DIR, "DE_FindAllMarkers_lane1_wilcox.csv")
A375_ETA_CSV <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_A375_flexv2.csv")
A375_DE_CSV  <- file.path(RESULTS_DIR, "DE_A375_flexv2_FindAllMarkers_wilcox.csv")

for (f in c(VCC_ETA_CSV, VCC_DE_CSV, A375_ETA_CSV, A375_DE_CSV)) {
  if (!file.exists(f)) stop("Required file not found: ", f,
                             "\nRun compute scripts first.")
}
if (!file.exists(VCC_MANIFEST_CSV)) {
  stop("VCC probe manifest not found: ", VCC_MANIFEST_CSV,
       "\nDownload and place at: ", VCC_MANIFEST_CSV,
       "\nSee data/README.md for download instructions.")
}
if (!file.exists(A375_MANIFEST_CSV)) {
  stop("A375 v2 probe manifest not found: ", A375_MANIFEST_CSV,
       "\nDownload and place at: ", A375_MANIFEST_CSV,
       "\nSee data/README.md for download instructions.")
}

# ==============================================================================
# VCC (Flex v1) probe counts
# ==============================================================================

message("Processing VCC (Flex v1) data...")

manifest_vcc <- read.csv(VCC_MANIFEST_CSV, comment.char = "#")
manifest_vcc <- manifest_vcc[manifest_vcc$included == TRUE, ]

# Parse gene name from probe_id (format: ENSG...|GENENAME|hash)
manifest_vcc$gene_name <- sapply(strsplit(manifest_vcc$probe_id, "\\|"), `[`, 2)

probe_counts_vcc <- manifest_vcc %>%
  group_by(gene_name) %>%
  summarise(n_probes = n_distinct(probe_id), .groups = "drop")

# Expressed genes (>= 10% of NTC cells)
vcc_eta <- read.csv(VCC_ETA_CSV)
expressed_vcc <- vcc_eta$gene[vcc_eta$pct_expressed >= MIN_PCT_EXPR]

probe_counts_vcc <- probe_counts_vcc[
  probe_counts_vcc$gene_name %in% expressed_vcc, ]

# Wilcox DE genes (any probe barcode, Lane 1)
vcc_de <- read.csv(VCC_DE_CSV)
de_genes_vcc <- unique(vcc_de$gene[
  !is.na(vcc_de$p_val_adj) & vcc_de$p_val_adj < PVAL_CUTOFF &
  abs(vcc_de$avg_log2FC) >= LOG2FC_CUTOFF
])

probe_counts_vcc$is_de <- ifelse(
  probe_counts_vcc$gene_name %in% de_genes_vcc, "DE", "Not DE")

message(sprintf("VCC: %d expressed genes  |  DE: %d  |  Not DE: %d",
  nrow(probe_counts_vcc),
  sum(probe_counts_vcc$is_de == "DE"),
  sum(probe_counts_vcc$is_de == "Not DE")))

# ==============================================================================
# A375 (Flex v2) probe counts
# ==============================================================================

message("Processing A375 (Flex v2) data...")

manifest_a375 <- read.csv(A375_MANIFEST_CSV, comment.char = "#")
manifest_a375 <- manifest_a375[manifest_a375$included == TRUE, ]

# Flex v2 manifest includes gene_name column directly
probe_counts_a375 <- manifest_a375 %>%
  group_by(gene_name) %>%
  summarise(n_probes = n_distinct(probe_id), .groups = "drop")

# Expressed genes (>= 10% of A375 cells)
a375_eta <- read.csv(A375_ETA_CSV)
expressed_a375 <- a375_eta$gene[a375_eta$pct_expressed >= MIN_PCT_EXPR]

probe_counts_a375 <- probe_counts_a375[
  probe_counts_a375$gene_name %in% expressed_a375, ]

# Wilcox DE genes from A375 FindAllMarkers (all cells, by ligand condition)
# NOTE: The published figure used PBS-only cells; see header for details.
a375_de <- read.csv(A375_DE_CSV)
de_genes_a375 <- unique(a375_de$gene[
  !is.na(a375_de$p_val_adj) & a375_de$p_val_adj < PVAL_CUTOFF &
  abs(a375_de$avg_log2FC) >= LOG2FC_CUTOFF
])

probe_counts_a375$is_de <- ifelse(
  probe_counts_a375$gene_name %in% de_genes_a375, "DE", "Not DE")

message(sprintf("A375: %d expressed genes  |  DE: %d  |  Not DE: %d",
  nrow(probe_counts_a375),
  sum(probe_counts_a375$is_de == "DE"),
  sum(probe_counts_a375$is_de == "Not DE")))

# ==============================================================================
# Plot side-by-side histograms
# ==============================================================================

make_panel <- function(probe_counts, panel_title) {
  probe_counts$is_de <- factor(probe_counts$is_de, levels = c("Not DE", "DE"))
  ggplot(probe_counts, aes(x = n_probes, fill = is_de)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1) +
    facet_wrap(~ is_de, ncol = 1) +
    scale_fill_manual(values = c("Not DE" = "grey", "DE" = "tomato3"),
                      name = NULL) +
    labs(x = "Number of probes per gene", y = "Density",
         title = panel_title) +
    theme_classic() +
    theme(
      plot.title           = element_text(size = 17, hjust = 0.5, color = "black"),
      axis.title           = element_text(size = 17, color = "black"),
      axis.text            = element_text(size = 14, color = "black"),
      axis.ticks           = element_line(color = "black"),
      axis.ticks.length    = unit(0.15, "cm"),
      axis.line            = element_line(color = "black", linewidth = 0.6),
      panel.border         = element_blank(),
      panel.grid.major     = element_blank(),
      panel.grid.minor     = element_blank(),
      strip.text           = element_blank(),
      strip.background     = element_blank(),
      legend.position      = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.text          = element_text(size = 14, color = "black")
    )
}

p_vcc  <- make_panel(probe_counts_vcc,  "Flex v1")
p_a375 <- make_panel(probe_counts_a375, "Flex v2")

p_combined <- p_vcc + p_a375 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "fig2d_probe_number_histogram_v1_vs_v2.pdf"),
       p_combined, width = 10, height = 5, device = cairo_pdf)
message("Saved: fig2d_probe_number_histogram_v1_vs_v2.pdf")

message("Done.")
print(sessionInfo())
