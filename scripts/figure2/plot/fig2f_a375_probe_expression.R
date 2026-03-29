#!/usr/bin/env Rscript
# ==============================================================================
# fig2f_a375_probe_expression.R
# Flex v1 DE genes — per-probe expression by probe barcode in A375 FLEXv2 (Fig 2F)
# ==============================================================================
#
# Takes the Flex v1 DE gene list (genes with eta² >= 1% in the VCC dataset,
# Figure 1F) and plots per-probe normalised expression distributions across
# probe barcodes in the A375 FLEXv2 PBS-cell dataset. This tests whether the
# same genes that show probe barcode effects in Flex v1 also show them in
# Flex v2, and whether individual probes drive the effect in each kit version.
#
# The paper features SRP9 and SPCS1. By default this script plots only those
# two genes; set GENES_TO_PLOT <- NULL to plot all Flex v1 DE genes that have
# >= 2 probes in the Flex v2 probe set.
#
# Note: The Flex v2 PBS-cell object uses `probe_barcode_number` (e.g. "A-A01")
# as the grouping variable, not the ligand condition name.
#
# INPUT
#   results/a375_probe_expression_cache.rds           — from 12_compute_a375_probe_expression.R
#   results/variance_explained_probe_barcode_NTC.csv  — from 04_compute_variance_explained.R
#                                                         (provides Flex v1 DE gene list)
#
# OUTPUT
#   results/figures/fig2f_SRP9_probe_expression_a375.pdf
#   results/figures/fig2f_SPCS1_probe_expression_a375.pdf
#   (and one PDF per additional gene if GENES_TO_PLOT is expanded)
#
# USAGE
#   Rscript scripts/figure2/plot/fig2f_a375_probe_expression.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

CACHE_RDS  <- file.path(RESULTS_DIR, "a375_probe_expression_cache.rds")
ETA_CSV    <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_NTC.csv")
ETA_THRESH <- 0.01   # Flex v1 DE threshold: eta² >= 1%
MIN_PROBES <- 2      # minimum probes per gene in Flex v2 to plot

# Set to NULL to plot all Flex v1 DE genes representable in Flex v2
GENES_TO_PLOT <- c("SRP9", "SPCS1")

for (f in c(CACHE_RDS, ETA_CSV)) {
  if (!file.exists(f)) stop("Required file not found: ", f,
                             "\nRun compute scripts first.")
}

# ==============================================================================
# Load data
# ==============================================================================

message("Loading A375 probe expression cache (PBS cells)...")
cache       <- readRDS(CACHE_RDS)
expr_data   <- cache$data      # normalised, probes x cells (Flex v2 probe IDs)
barcodes    <- cache$barcodes  # probe_barcode_number per cell
all_probes  <- rownames(expr_data)
rm(cache); gc()

bc_levels <- sort(unique(barcodes))
message(sprintf("Probes: %d  |  Cells: %d  |  Probe barcodes: %d",
                nrow(expr_data), ncol(expr_data), length(bc_levels)))

# Flex v1 DE gene list (eta² >= 1% in VCC NTC cells)
v1_de_genes <- read.csv(ETA_CSV) %>%
  filter(eta_sq_bc >= ETA_THRESH) %>%
  pull(gene)
message(sprintf("Flex v1 DE genes loaded: %d", length(v1_de_genes)))

# ==============================================================================
# Build probe map in Flex v2 data for Flex v1 DE genes
# ==============================================================================

build_probe_map <- function(gene_list) {
  pm <- lapply(gene_list, function(gene) {
    grep(paste0("-", gene, "-"), all_probes, value = TRUE)
  })
  names(pm) <- gene_list
  Filter(function(x) length(x) >= MIN_PROBES, pm)
}

if (is.null(GENES_TO_PLOT)) {
  probe_map <- build_probe_map(v1_de_genes)
} else {
  probe_map <- build_probe_map(GENES_TO_PLOT)
  missing <- setdiff(GENES_TO_PLOT, names(probe_map))
  if (length(missing) > 0) {
    message("Genes not found / < ", MIN_PROBES,
            " probes in Flex v2: ", paste(missing, collapse = ", "))
  }
}
message(sprintf("Flex v1 DE genes with >= %d probes in Flex v2 data: %d",
                MIN_PROBES, length(probe_map)))

if (length(probe_map) == 0) {
  stop("No genes to plot. Check that GENES_TO_PLOT are in the Flex v1 DE list ",
       "and have >= ", MIN_PROBES, " probes in the Flex v2 probe set.")
}

# ==============================================================================
# Subset expression and build long-format data frame
# ==============================================================================

target_probes <- unique(unlist(probe_map))
expr_sub      <- as.matrix(expr_data[target_probes, ])
rm(expr_data); gc()

PROBE_COLORS <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

df_long <- as.data.frame(t(expr_sub), check.names = FALSE)
df_long$probe_barcode <- factor(barcodes, levels = bc_levels)
df_long <- pivot_longer(df_long, cols = -probe_barcode,
                         names_to = "probe_id", values_to = "expression") %>%
  mutate(
    probe_label = substr(probe_id, nchar(probe_id) - 6, nchar(probe_id)),
    gene        = sub(".*\\|([^|]+)\\|[^|]+$", "\\1", probe_id)
  )

# ==============================================================================
# Plot function
# ==============================================================================

plot_gene <- function(gene, probe_map, df_long, out_dir) {
  probes   <- probe_map[[gene]]
  n_probes <- length(probes)
  out_file <- file.path(out_dir, sprintf("fig2f_%s_probe_expression_a375.pdf", gene))

  df_gene  <- filter(df_long, gene == !!gene)
  colors   <- PROBE_COLORS[seq_len(n_probes)]

  p <- ggplot(df_gene, aes(x = probe_barcode, y = expression,
                             fill = probe_label)) +
    geom_violin(position = position_dodge(0.9),
                scale = "width", trim = TRUE) +
    scale_fill_manual(values = setNames(colors, unique(df_gene$probe_label)),
                      name = "Probe") +
    labs(title = bquote(italic(.(gene))),
         x = "Probe barcode", y = "Normalised expression") +
    theme_classic() +
    theme(
      plot.title   = element_text(size = 18, hjust = 0.5, color = "black"),
      axis.title   = element_text(size = 16, color = "black"),
      axis.text    = element_text(size = 14, color = "black"),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.ticks   = element_line(color = "black"),
      axis.line    = element_line(color = "black", linewidth = 0.6),
      legend.text  = element_text(size = 13, color = "black"),
      legend.title = element_text(size = 13, color = "black")
    )

  ggsave(out_file, plot = p, width = 10, height = 3, device = cairo_pdf)
  message("Saved: ", basename(out_file))
}

for (gene in names(probe_map)) {
  plot_gene(gene, probe_map, df_long, FIGURES_DIR)
}

message("Done.")
print(sessionInfo())
