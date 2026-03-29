#!/usr/bin/env Rscript
# ==============================================================================
# 05_compute_fdp_sweep.R
# False discovery proportion sweep — confounded probe barcode DE (Figures 1G, 1H)
# ==============================================================================
#
# Quantifies the false positive rate that arises when guide identity is
# confounded with probe barcode identity in a case-control comparison.
#
# For each CRISPR target gene in the VCC dataset, two DE analyses are run:
#
#   Ground truth  — NTC cells (TRUTH_BC) vs target gene cells (TRUTH_BC)
#     Both groups from the same probe barcode: no confounding.
#
#   Confounded    — target gene cells (CONFOUNDED_BC) vs NTC cells (TRUTH_BC)
#     Target and control come from different probe barcodes, simulating an
#     experiment where biological sample identity is confounded with barcode.
#
# For each gene, DE results are merged and FDP = FP / (FP + TP) is computed
# across a sweep of |log2FC| thresholds (0 to 2). Results are saved per gene
# and combined into a single metrics table for plotting.
#
# Previously computed per-gene CSVs are cached and skipped on re-runs.
#
# Barcode pair used (represents average technical effect among the 16 barcodes):
#   TRUTH_BC = BC005, CONFOUNDED_BC = BC008
#
# INPUT
#   results/seurat_vcc_annotated.qs  — from 01_setup_vcc_seurat.R
#
# OUTPUT
#   results/fdp_sweep_BC005_vs_BC008/
#     de_truth_<gene>.csv          — ground truth DE per target gene
#     de_confounded_<gene>.csv     — confounded DE per target gene
#     metrics_<gene>.csv           — FDP metrics per gene across |log2FC| sweep
#     metrics_all_genes.csv        — all genes stacked (used for Figure 1H)
#     gene_status.csv              — processing status for every target gene
#
# USAGE
#   Rscript scripts/figure1/compute/05_compute_fdp_sweep.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
})

TRUTH_BC      <- "BC005"
CONFOUNDED_BC <- "BC008"
NTC_LABEL     <- "non-targeting"
MIN_CELLS     <- 20
MIN_PCT       <- 0.1
LOG2FC_CUTOFFS <- seq(0, 2, by = 0.1)
PVAL_CUTOFF   <- 0.05

RESULTS_DIR  <- "results"
SWEEP_DIR    <- file.path(RESULTS_DIR, paste0("fdp_sweep_", TRUTH_BC, "_vs_", CONFOUNDED_BC))
VCC_QS       <- file.path(RESULTS_DIR, "seurat_vcc_annotated.qs")
COMBINED_CSV <- file.path(SWEEP_DIR, "metrics_all_genes.csv")

if (!file.exists(VCC_QS)) {
  stop("VCC Seurat not found: ", VCC_QS, "\nRun 01_setup_vcc_seurat.R first.")
}

dir.create(SWEEP_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Helpers
# ==============================================================================

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", gsub("[/\\\\]", "-", x))

calc_metrics <- function(lfc_truth, padj_truth, lfc_pred, padj_pred,
                         lfc_cutoff, pval_cutoff) {
  gt   <- abs(lfc_truth) >= lfc_cutoff & padj_truth < pval_cutoff
  pred <- abs(lfc_pred)  >= lfc_cutoff & padj_pred  < pval_cutoff
  TP   <- sum(gt &  pred, na.rm = TRUE)
  FP   <- sum(!gt & pred, na.rm = TRUE)
  FN   <- sum(gt & !pred, na.rm = TRUE)
  FDP  <- if ((TP + FP) > 0) FP / (TP + FP) else NA_real_
  data.frame(TP = TP, FP = FP, FN = FN, FDP = FDP)
}

# ==============================================================================
# Load Seurat
# ==============================================================================

message(sprintf("=== FDP sweep: %s (truth) vs %s (confounded) ===\n",
                TRUTH_BC, CONFOUNDED_BC))

message("Loading annotated VCC Seurat object...")
seurat_vcc <- qread(VCC_QS)
message(sprintf("Loaded: %d cells", ncol(seurat_vcc)))

cells_ntc <- WhichCells(seurat_vcc,
  expression = target_gene == NTC_LABEL & probe_barcode == TRUTH_BC)
message(sprintf("NTC cells in %s: %d", TRUTH_BC, length(cells_ntc)))

if (length(cells_ntc) < MIN_CELLS) {
  stop("Too few NTC cells in TRUTH_BC = ", TRUTH_BC)
}

# ==============================================================================
# Identify testable target genes (present in both barcodes)
# ==============================================================================

targets_truth <- unique(seurat_vcc$target_gene[
  seurat_vcc$probe_barcode == TRUTH_BC &
  seurat_vcc$target_gene   != NTC_LABEL &
  !is.na(seurat_vcc$target_gene)])

targets_conf  <- unique(seurat_vcc$target_gene[
  seurat_vcc$probe_barcode == CONFOUNDED_BC &
  seurat_vcc$target_gene   != NTC_LABEL &
  !is.na(seurat_vcc$target_gene)])

targets <- sort(intersect(targets_truth, targets_conf))
message(sprintf("Target genes testable in both barcodes: %d\n", length(targets)))

# ==============================================================================
# Per-gene FDP sweep
# ==============================================================================

status_df    <- data.frame(target_gene = targets, status = NA_character_,
                           stringsAsFactors = FALSE)
results_list <- list()

for (i in seq_along(targets)) {
  tg      <- targets[i]
  tg_safe <- safe_name(tg)
  message(sprintf("[%d/%d] %s", i, length(targets), tg))

  truth_csv <- file.path(SWEEP_DIR, paste0("de_truth_",      tg_safe, ".csv"))
  conf_csv  <- file.path(SWEEP_DIR, paste0("de_confounded_", tg_safe, ".csv"))
  met_csv   <- file.path(SWEEP_DIR, paste0("metrics_",       tg_safe, ".csv"))

  cells_truth <- WhichCells(seurat_vcc,
    expression = target_gene == tg & probe_barcode == TRUTH_BC)
  cells_conf  <- WhichCells(seurat_vcc,
    expression = target_gene == tg & probe_barcode == CONFOUNDED_BC)

  if (length(cells_truth) < MIN_CELLS || length(cells_conf) < MIN_CELLS) {
    status_df$status[i] <- "skipped_few_cells"
    next
  }

  # Load cached metrics if available
  if (file.exists(met_csv)) {
    results_list[[tg]] <- read.csv(met_csv)
    status_df$status[i] <- "cached"
    next
  }

  # Ground truth DE
  if (!file.exists(truth_csv)) {
    de_truth <- tryCatch(
      FindMarkers(seurat_vcc, ident.1 = cells_ntc, ident.2 = cells_truth,
                  group.by = NULL, logfc.threshold = 0, min.pct = MIN_PCT,
                  test.use = "wilcox"),
      error = function(e) NULL
    )
    if (is.null(de_truth)) { status_df$status[i] <- "error"; next }
    de_truth$gene <- rownames(de_truth)
    de_truth$p_val_adj[is.na(de_truth$p_val_adj)] <- 1
    write.csv(de_truth, truth_csv, row.names = FALSE)
  } else {
    de_truth <- read.csv(truth_csv)
  }

  # Confounded DE
  if (!file.exists(conf_csv)) {
    de_conf <- tryCatch(
      FindMarkers(seurat_vcc, ident.1 = cells_conf, ident.2 = cells_ntc,
                  group.by = NULL, logfc.threshold = 0, min.pct = MIN_PCT,
                  test.use = "wilcox"),
      error = function(e) NULL
    )
    if (is.null(de_conf)) { status_df$status[i] <- "error"; next }
    de_conf$gene <- rownames(de_conf)
    de_conf$p_val_adj[is.na(de_conf$p_val_adj)] <- 1
    write.csv(de_conf, conf_csv, row.names = FALSE)
  } else {
    de_conf <- read.csv(conf_csv)
  }

  # Merge and sweep
  df_merged <- merge(
    de_truth[, c("gene", "avg_log2FC", "p_val_adj")],
    de_conf[,  c("gene", "avg_log2FC", "p_val_adj")],
    by = "gene", suffixes = c("_truth", "_confounded")
  )
  if (nrow(df_merged) == 0) { status_df$status[i] <- "skipped_no_overlap"; next }

  df_merged$p_val_adj_truth[is.na(df_merged$p_val_adj_truth)]           <- 1
  df_merged$p_val_adj_confounded[is.na(df_merged$p_val_adj_confounded)] <- 1

  metrics_rows <- lapply(LOG2FC_CUTOFFS, function(lfc_cut) {
    m <- calc_metrics(df_merged$avg_log2FC_truth,      df_merged$p_val_adj_truth,
                      df_merged$avg_log2FC_confounded, df_merged$p_val_adj_confounded,
                      lfc_cut, PVAL_CUTOFF)
    cbind(data.frame(target_gene   = tg,
                     log2fc_cutoff = lfc_cut,
                     pval_cutoff   = PVAL_CUTOFF,
                     n_genes_total = nrow(df_merged),
                     n_cells_truth = length(cells_truth),
                     n_cells_conf  = length(cells_conf),
                     n_cells_ntc   = length(cells_ntc)),
          m)
  })

  metrics_df <- do.call(rbind, metrics_rows)
  write.csv(metrics_df, met_csv, row.names = FALSE)
  results_list[[tg]] <- metrics_df
  status_df$status[i] <- "completed"
}

# ==============================================================================
# Combine and save
# ==============================================================================

write.csv(status_df, file.path(SWEEP_DIR, "gene_status.csv"), row.names = FALSE)
message("\nStatus summary:")
print(table(status_df$status))

if (length(results_list) == 0) stop("No genes were successfully processed.")

df_all <- do.call(rbind, results_list)
rownames(df_all) <- NULL
write.csv(df_all, COMBINED_CSV, row.names = FALSE)
message(sprintf("\nSaved combined metrics (%d rows) to:\n  %s", nrow(df_all), COMBINED_CSV))

message("Done.")
print(sessionInfo())
