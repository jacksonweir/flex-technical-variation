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
# The paper features SRP9 and SPCS1. Set GENES_TO_PLOT <- NULL to plot all
# Flex v1 DE genes that have >= MIN_PROBES probes in the Flex v2 probe set.
#
# NOTE: Requires Matrix >= 1.7 (R 4.5+). Older Matrix versions may drop
# dgCMatrix rownames when loading the probe expression cache.
#
# INPUT
#   results/a375_probe_expression_cache.rds           — from 12_compute_a375_probe_expression.R
#   results/variance_explained_probe_barcode_NTC.csv  — from 04_compute_variance_explained.R
#                                                         (provides Flex v1 DE gene list)
#
#   The cache corresponds to flex_v2_probe_expr_sparse_cache_10k_feat_filt.rds
#   from the Zenodo deposit (https://doi.org/10.5281/zenodo.19363777). Download
#   and process with 12_compute_a375_probe_expression.R, which
#   saves a copy to results/a375_probe_expression_cache.rds.
#
# OUTPUT
#   results/figures/fig2f_SRP9_probe_expression_a375.pdf
#   results/figures/fig2f_SPCS1_probe_expression_a375.pdf
#
# USAGE
#   Rscript scripts/figure2/plot/fig2f_a375_probe_expression.R
#
# PACKAGES
#   ggplot2, dplyr, tidyr, Matrix, dittoSeq, scales
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(dittoSeq)
  library(scales)
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

PROBE_COLORS <- dittoColors()

for (f in c(CACHE_RDS, ETA_CSV)) {
  if (!file.exists(f)) stop("Required file not found: ", f,
                             "\nRun compute scripts first.")
}

# ==============================================================================
# Load data
# ==============================================================================

message("Loading A375 Flex v2 probe expression cache (PBS cells)...")
cache       <- readRDS(CACHE_RDS)
expr_data   <- cache$data      # normalised, probes x cells (Flex v2 probe IDs)
barcodes    <- cache$barcodes  # probe_barcode_number per cell
all_probes  <- rownames(expr_data)
rm(cache); gc()

bc_levels <- sort(unique(barcodes))
message(sprintf("Probes: %d  |  Cells: %d  |  Probe barcodes: %d",
                nrow(expr_data), ncol(expr_data), length(bc_levels)))

# Flex v1 DE gene list (eta² >= 1% in VCC NTC cells)
v1_de_genes <- read.csv(ETA_CSV) |>
  filter(eta_sq_bc >= ETA_THRESH) |>
  pull(gene)
message(sprintf("Flex v1 DE genes loaded: %d", length(v1_de_genes)))

# ==============================================================================
# Build probe map in Flex v2 data for Flex v1 DE genes
# Probe ID format: ENSG...-GENENAME-hash (dash separator in these caches)
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
for (g in names(probe_map)) message("  ", g, ": ", paste(probe_map[[g]], collapse = ", "))

if (length(probe_map) == 0) {
  stop("No genes to plot. Check GENES_TO_PLOT are in the Flex v1 DE list ",
       "and have >= ", MIN_PROBES, " probes in the Flex v2 probe set.")
}

# ==============================================================================
# Subset expression and build long-format data frame
# ==============================================================================

target_probes <- unique(unlist(probe_map))
expr_sub      <- as.matrix(expr_data[target_probes, ])
rm(expr_data); gc()

df_long <- as.data.frame(t(expr_sub), check.names = FALSE)
df_long$probe_barcode <- factor(barcodes, levels = bc_levels)
df_long <- pivot_longer(df_long, cols = -probe_barcode,
                         names_to = "probe_id", values_to = "expression") |>
  mutate(
    probe_label = substr(probe_id, nchar(probe_id) - 6, nchar(probe_id)),
    gene        = sub("^ENSG[0-9]+-([^-]+)-.*", "\\1", probe_id)
  )

# ==============================================================================
# Plot
# ==============================================================================

for (gene in names(probe_map)) {
  n_probes <- length(probe_map[[gene]])
  df_gene  <- filter(df_long, gene == !!gene)
  colors   <- setNames(PROBE_COLORS[seq_len(n_probes)], unique(df_gene$probe_label))
  out_file <- file.path(FIGURES_DIR, sprintf("fig2f_%s_probe_expression_a375.pdf", gene))

  p <- ggplot(df_gene, aes(x = probe_barcode, y = expression, fill = probe_label)) +
    geom_violin(position = position_dodge(0.9), scale = "width",
                trim = TRUE, linewidth = 0.3) +
    scale_fill_manual(values = colors, name = "Probe ID") +
    scale_y_continuous(breaks = breaks_width(1), minor_breaks = NULL) +
    labs(title = bquote(italic(.(gene))),
         x = "Probe barcode", y = "Expression") +
    theme_classic() +
    theme(
      text          = element_text(color = "black"),
      plot.title    = element_text(size = 26, hjust = 0.5, color = "black",
                                   face = "italic"),
      axis.title    = element_text(size = 24, color = "black"),
      axis.text     = element_text(size = 20, color = "black"),
      axis.text.x   = element_text(size = 20, angle = 45, hjust = 1),
      axis.ticks    = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line     = element_line(color = "black", linewidth = 0.6),
      panel.border  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text   = element_text(size = 20, color = "black"),
      legend.title  = element_text(size = 22, color = "black")
    )

  ggsave(out_file, plot = p, width = 10, height = 3, device = cairo_pdf)
  message("Saved: ", basename(out_file))
}

message("Done.")
print(sessionInfo())
