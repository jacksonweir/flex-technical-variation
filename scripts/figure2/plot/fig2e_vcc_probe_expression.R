#!/usr/bin/env Rscript
# ==============================================================================
# fig2e_vcc_probe_expression.R
# Per-probe expression by probe barcode — VCC Flex v1 DE genes (Figure 2E)
# ==============================================================================
#
# For each Flex v1 DE gene with at least two probes (e.g. SRP9, SPCS1), plots
# violin distributions of normalised expression for each probe separately,
# grouped by probe barcode. This reveals whether the overall gene-level
# barcode effect is driven by one probe or distributed across all probes
# targeting the gene.
#
# The paper features SRP9 (highest eta² gene) and SPCS1 (second highest).
# Set GENES_TO_PLOT <- NULL to plot all DE genes with >= MIN_PROBES probes.
#
# NOTE: Requires Matrix >= 1.7 (R 4.5+). Older Matrix versions may drop
# dgCMatrix rownames when loading the probe expression cache.
#
# INPUT
#   results/vcc_probe_expression_cache.rds        — from 11_compute_vcc_probe_expression.R
#   results/variance_explained_probe_barcode_NTC.csv  — from 04_compute_variance_explained.R
#
#   The cache file corresponds to vcc_probe_expr_sparse_cache_10k_feat_filt.rds
#   from the Zenodo deposit (https://doi.org/10.5281/zenodo.19363777). Download
#   and process with 11_compute_vcc_probe_expression.R, which
#   saves a copy to results/vcc_probe_expression_cache.rds.
#
# OUTPUT
#   results/figures/fig2e_SRP9_probe_expression_vcc.pdf
#   results/figures/fig2e_SPCS1_probe_expression_vcc.pdf
#
# USAGE
#   Rscript scripts/figure2/plot/fig2e_vcc_probe_expression.R
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

CACHE_RDS  <- file.path(RESULTS_DIR, "vcc_probe_expression_cache.rds")
ETA_CSV    <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_NTC.csv")
ETA_THRESH <- 0.01   # eta² >= 1% defines a DE gene
MIN_PROBES <- 2      # minimum probes per gene to plot

# Set to NULL to plot all DE genes with >= MIN_PROBES probes
GENES_TO_PLOT <- c("SRP9", "SPCS1")

PROBE_COLORS <- dittoColors()

for (f in c(CACHE_RDS, ETA_CSV)) {
  if (!file.exists(f)) stop("Required file not found: ", f,
                             "\nRun compute scripts first.")
}

# ==============================================================================
# Load data
# ==============================================================================

message("Loading VCC probe expression cache...")
cache       <- readRDS(CACHE_RDS)
expr_data   <- cache$data      # normalised, probes x cells
barcodes    <- cache$barcodes
all_probes  <- rownames(expr_data)
rm(cache); gc()

bc_levels <- paste0("BC", sprintf("%03d", 1:16))
bc_levels <- bc_levels[bc_levels %in% unique(barcodes)]
message(sprintf("Probes: %d  |  Cells: %d  |  Probe barcodes: %d",
                nrow(expr_data), ncol(expr_data), length(bc_levels)))

# DE gene list from eta-squared results
de_genes <- read.csv(ETA_CSV) |>
  filter(eta_sq_bc >= ETA_THRESH) |>
  pull(gene)
message(sprintf("DE genes (eta\u00b2 >= 1%%): %d", length(de_genes)))

# ==============================================================================
# Build probe map: gene → probes targeting it (requires >= MIN_PROBES)
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
  probe_map <- build_probe_map(de_genes)
} else {
  genes_use <- intersect(GENES_TO_PLOT, de_genes)
  if (length(genes_use) == 0) stop("None of GENES_TO_PLOT found in DE gene list.")
  if (length(genes_use) < length(GENES_TO_PLOT)) {
    message("Note: not in DE list: ",
            paste(setdiff(GENES_TO_PLOT, genes_use), collapse = ", "))
  }
  probe_map <- build_probe_map(genes_use)
}
message(sprintf("Genes with >= %d probes to plot: %d", MIN_PROBES, length(probe_map)))
for (g in names(probe_map)) message("  ", g, ": ", paste(probe_map[[g]], collapse = ", "))

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
# Plot function
# ==============================================================================

for (gene in names(probe_map)) {
  n_probes <- length(probe_map[[gene]])
  df_gene  <- filter(df_long, gene == !!gene)
  colors   <- setNames(PROBE_COLORS[seq_len(n_probes)], unique(df_gene$probe_label))
  out_file <- file.path(FIGURES_DIR, sprintf("fig2e_%s_probe_expression_vcc.pdf", gene))

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
