#!/usr/bin/env Rscript
# ==============================================================================
# 10_compute_variance_explained_a375.R
# Eta-squared (variance explained) by probe barcode — A375 FLEXv2 (Figure 2D)
# ==============================================================================
#
# Computes the fraction of total expression variance attributable to probe
# barcode identity (= ligand treatment condition) for each expressed gene in
# the A375 FLEXv2 dataset: eta² = SS_between / SS_total from a one-way ANOVA.
#
# Unlike the VCC analysis (script 04), all cells are included here — there is
# no NTC-only subset because each probe barcode corresponds to a distinct
# biological condition. Consequently, eta² reflects both biological (ligand)
# and technical (barcode) effects.
#
# Gene filter: expressed in >= 5% of cells.
# SS computed via sparse matrix algebra — no matrix densification required.
# P-values from the F-distribution, FDR-corrected (BH).
#
# INPUT
#   results/seurat_a375_flexv2_filtered.qs  — from 08_setup_a375_seurat.R
#
# OUTPUT
#   results/variance_explained_probe_barcode_A375_flexv2.csv
#     Columns: gene, eta_sq_bc, F_stat, p_val, p_adj, pct_expressed, rank_bc
#
# USAGE
#   Rscript scripts/figure2/compute/10_compute_variance_explained_a375.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(Matrix)
  library(dplyr)
  library(readr)
})

RESULTS_DIR  <- "results"
IN_QS        <- file.path(RESULTS_DIR, "seurat_a375_flexv2_filtered.qs")
OUT_CSV      <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_A375_flexv2.csv")
MIN_PCT_EXPR <- 0.05

if (!file.exists(IN_QS)) {
  stop("Filtered A375 Seurat not found: ", IN_QS,
       "\nRun 08_setup_a375_seurat.R first.")
}

if (file.exists(OUT_CSV)) {
  message("Already exists, skipping: ", OUT_CSV)
  quit(status = 0)
}

# ==============================================================================
# Load and extract expression matrix
# ==============================================================================

message("Loading filtered A375 FLEXv2 Seurat object...")
obj <- qread(IN_QS)
message(sprintf("Cells: %d", ncol(obj)))

bc_vec <- obj$probe_barcode
expr   <- GetAssayData(obj, layer = "data")   # genes x cells, sparse
message(sprintf("Matrix: %d genes x %d cells", nrow(expr), ncol(expr)))
rm(obj); gc()

pct_expr   <- Matrix::rowMeans(expr > 0)
keep_genes <- pct_expr >= MIN_PCT_EXPR
expr_filt  <- expr[keep_genes, ]
message(sprintf("Genes kept (>= %.0f%% cells expressing): %d / %d",
                MIN_PCT_EXPR * 100, nrow(expr_filt), nrow(expr)))
rm(expr); gc()

# ==============================================================================
# Eta-squared via sparse matrix algebra
# ==============================================================================

message("Computing eta-squared for probe barcode...")
groups <- sort(unique(bc_vec))
K      <- length(groups)
N      <- length(bc_vec)
df1    <- K - 1
df2    <- N - K

grand_mean <- Matrix::rowMeans(expr_filt)
SS_total   <- as.numeric(Matrix::rowSums(expr_filt^2) - N * grand_mean^2)
SS_between <- rep(0.0, nrow(expr_filt))

for (g in groups) {
  idx <- which(bc_vec == g)
  if (length(idx) < 2) next
  gm         <- as.numeric(Matrix::rowMeans(expr_filt[, idx, drop = FALSE]))
  SS_between <- SS_between + length(idx) * gm^2
}
SS_between <- SS_between - N * as.numeric(grand_mean)^2

eta_sq           <- SS_between / SS_total
eta_sq[SS_total < 1e-10] <- 0
eta_sq[eta_sq < 0]       <- 0
eta_sq[eta_sq > 1]       <- 1

SS_within <- pmax(SS_total - SS_between, 0)
F_stat    <- (SS_between / df1) / (SS_within / df2)
F_stat[SS_total < 1e-10] <- 0

p_val <- pf(F_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
p_adj <- p.adjust(p_val, method = "BH")

# ==============================================================================
# Assemble and save results
# ==============================================================================

df_eta <- data.frame(
  gene          = rownames(expr_filt),
  eta_sq_bc     = eta_sq,
  F_stat        = F_stat,
  p_val         = p_val,
  p_adj         = p_adj,
  pct_expressed = as.numeric(pct_expr[keep_genes]),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(eta_sq_bc)) %>%
  mutate(rank_bc = seq_len(n()))

message(sprintf("\nGenes with eta²_bc >= 1%%:   %d", sum(df_eta$eta_sq_bc >= 0.01)))
message(sprintf("Genes with eta²_bc >= 5%%:   %d", sum(df_eta$eta_sq_bc >= 0.05)))
message("Top 10 genes by eta²:")
print(head(df_eta[, c("gene", "eta_sq_bc", "rank_bc")], 10))

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(df_eta, OUT_CSV)
message("Saved: ", OUT_CSV)

message("Done.")
print(sessionInfo())
