#!/usr/bin/env Rscript
# ==============================================================================
# 04_compute_variance_explained.R
# Eta-squared (variance explained) by probe barcode and lane — NTC cells only
# ==============================================================================
#
# For each expressed gene, computes the fraction of total expression variance
# attributable to probe barcode identity (eta² = SS_between / SS_total from a
# one-way ANOVA). The same calculation is repeated with lane as the grouping
# factor for comparison.
#
# Only non-targeting control (NTC) cells are used. Because NTC cells carry no
# guide RNA perturbation, all cells are biologically equivalent by design. Any
# variance explained by probe barcode is therefore a purely technical artifact.
#
# Gene filter: expressed in >= 5% of NTC cells.
# SS computed via sparse matrix algebra — no densification required.
# P-values derived from the F-distribution and FDR-corrected (BH).
#
# INPUT
#   results/seurat_vcc_annotated.qs  — from 01_setup_vcc_seurat.R
#
# OUTPUT
#   results/variance_explained_probe_barcode_NTC.csv
#     Columns: gene, eta_sq_bc, F_stat_bc, p_val_bc, p_adj_bc,
#              eta_sq_lane, F_stat_lane, p_val_lane, p_adj_lane,
#              pct_expressed, rank_bc, rank_lane
#
# USAGE
#   Rscript scripts/figure1/compute/04_compute_variance_explained.R
#
# PACKAGES
#   Seurat, qs, Matrix, dplyr, readr
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(Matrix)
  library(dplyr)
  library(readr)
})

RESULTS_DIR  <- "results"
VCC_QS       <- file.path(RESULTS_DIR, "seurat_vcc_annotated.qs")
OUT_CSV      <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_NTC.csv")
NTC_LABEL    <- "non-targeting"
MIN_PCT_EXPR <- 0.05

if (!file.exists(VCC_QS)) {
  stop("VCC Seurat not found: ", VCC_QS, "\nRun 01_setup_vcc_seurat.R first.")
}

if (file.exists(OUT_CSV)) {
  message("Already exists, skipping: ", OUT_CSV)
  quit(status = 0)
}

# ==============================================================================
# Load and subset to NTC cells
# ==============================================================================

message("Loading annotated VCC Seurat object...")
seurat_vcc <- qread(VCC_QS)

cells_ntc <- colnames(seurat_vcc)[
  !is.na(seurat_vcc$target_gene) &
  seurat_vcc$target_gene == NTC_LABEL &
  !is.na(seurat_vcc$probe_barcode)
]
message(sprintf("NTC cells: %d / %d total", length(cells_ntc), ncol(seurat_vcc)))

seurat_ntc <- subset(seurat_vcc, cells = cells_ntc)
rm(seurat_vcc); gc()

message("NTC cells per probe barcode:")
print(sort(table(seurat_ntc$probe_barcode)))
message("NTC cells per lane:")
print(sort(table(seurat_ntc$lane)))

# ==============================================================================
# Extract expression matrix and filter genes
# ==============================================================================

message("Extracting log-normalised expression matrix...")
expr <- GetAssayData(seurat_ntc, layer = "data")   # genes x cells, sparse

pct_expr   <- Matrix::rowMeans(expr > 0)
keep_genes <- pct_expr >= MIN_PCT_EXPR
expr_filt  <- expr[keep_genes, ]
message(sprintf("Genes kept (>= %.0f%% cells expressing): %d / %d",
                MIN_PCT_EXPR * 100, nrow(expr_filt), nrow(expr)))

# ==============================================================================
# Eta-squared function (sparse matrix algebra)
# ==============================================================================

compute_eta_sq <- function(mat, group_vec) {
  # mat:       genes x cells sparse Matrix
  # group_vec: factor/character vector of length ncol(mat)
  # Returns:   data.frame with eta_sq, F_stat, p_val, p_adj per gene

  N      <- ncol(mat)
  groups <- sort(unique(group_vec))
  K      <- length(groups)
  df1    <- K - 1
  df2    <- N - K

  grand_mean <- Matrix::rowMeans(mat)

  SS_total   <- as.numeric(Matrix::rowSums(mat^2) - N * grand_mean^2)
  SS_between <- rep(0.0, nrow(mat))

  for (g in groups) {
    idx <- which(group_vec == g)
    if (length(idx) < 2) next
    gm          <- as.numeric(Matrix::rowMeans(mat[, idx, drop = FALSE]))
    SS_between  <- SS_between + length(idx) * gm^2
  }
  SS_between <- SS_between - N * as.numeric(grand_mean)^2

  eta_sq           <- SS_between / SS_total
  eta_sq[SS_total < 1e-10] <- 0
  eta_sq[eta_sq < 0]       <- 0
  eta_sq[eta_sq > 1]       <- 1

  SS_within        <- pmax(SS_total - SS_between, 0)
  F_stat           <- (SS_between / df1) / (SS_within / df2)
  F_stat[SS_total < 1e-10] <- 0

  p_val <- pf(F_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
  p_adj <- p.adjust(p_val, method = "BH")

  data.frame(gene = rownames(mat), eta_sq = eta_sq,
             F_stat = F_stat, p_val = p_val, p_adj = p_adj,
             stringsAsFactors = FALSE)
}

# ==============================================================================
# Compute eta-squared for probe barcode and lane
# ==============================================================================

message("Computing eta-squared for probe barcode (K=16)...")
res_bc   <- compute_eta_sq(expr_filt, seurat_ntc$probe_barcode)

message("Computing eta-squared for lane (K=3)...")
res_lane <- compute_eta_sq(expr_filt, seurat_ntc$lane)

# ==============================================================================
# Assemble results table and save
# ==============================================================================

df_eta <- data.frame(
  gene          = res_bc$gene,
  eta_sq_bc     = res_bc$eta_sq,
  F_stat_bc     = res_bc$F_stat,
  p_val_bc      = res_bc$p_val,
  p_adj_bc      = res_bc$p_adj,
  eta_sq_lane   = res_lane$eta_sq,
  F_stat_lane   = res_lane$F_stat,
  p_val_lane    = res_lane$p_val,
  p_adj_lane    = res_lane$p_adj,
  pct_expressed = as.numeric(pct_expr[keep_genes]),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(eta_sq_bc)) %>%
  mutate(rank_bc = seq_len(n())) %>%
  arrange(desc(eta_sq_lane)) %>%
  mutate(rank_lane = seq_len(n())) %>%
  arrange(rank_bc)

message(sprintf("\nGenes with eta²_bc >= 1%%:   %d", sum(df_eta$eta_sq_bc   >= 0.01)))
message(sprintf("Genes with eta²_lane >= 1%%: %d", sum(df_eta$eta_sq_lane >= 0.01)))
message("Top 5 genes by probe barcode eta²:")
print(head(df_eta[, c("gene", "eta_sq_bc", "rank_bc")], 5))

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(df_eta, OUT_CSV)
message("Saved: ", OUT_CSV)

message("Done.")
print(sessionInfo())
