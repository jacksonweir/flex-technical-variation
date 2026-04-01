#!/usr/bin/env Rscript
# ==============================================================================
# 00_prep_violin_data.R
# Extract per-cell expression of the top two eta²-ranked genes (SRP9, SPCS1)
# from NTC cells across all 3 VCC lanes for Figure 1D violin plots.
#
# INPUT
#   results/seurat_vcc_annotated.qs        — from 01_setup_vcc_seurat.R
#   results/variance_explained_probe_barcode_NTC.csv  — from 04_compute_variance_explained.R
#
# OUTPUT
#   results/vcc_ntc_violin_data.csv
#     long-format: cell, gene, expression, probe_barcode, lane
#
# NOTE: Run 01_setup_vcc_seurat.R and 04_compute_variance_explained.R first.
#       04 must complete before this script to confirm gene identities, but
#       you can also set TOP_GENES explicitly below if needed.
#
# USAGE
#   Rscript scripts/figure1/00_prep_violin_data.R
#
# PACKAGES
#   Seurat, qs, Matrix, dplyr, tidyr, tibble
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

RESULTS_DIR <- "results"
VCC_QS      <- file.path(RESULTS_DIR, "seurat_vcc_annotated.qs")
ETA_CSV     <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_NTC.csv")
OUT_CSV     <- file.path(RESULTS_DIR, "vcc_ntc_violin_data.csv")
NTC_LABEL   <- "non-targeting"
N_TOP       <- 2   # number of top eta² genes to extract

if (!file.exists(VCC_QS)) {
  stop("VCC annotated Seurat not found: ", VCC_QS,
       "\nRun 01_setup_vcc_seurat.R first.")
}
if (!file.exists(ETA_CSV)) {
  stop("Eta-squared CSV not found: ", ETA_CSV,
       "\nRun 04_compute_variance_explained.R first.")
}

# ==============================================================================
# Identify top N genes by probe barcode eta²
# ==============================================================================

eta      <- read.csv(ETA_CSV) |> dplyr::arrange(rank_bc)
top_genes <- eta$gene[seq_len(N_TOP)]
message(sprintf("Top %d eta² genes: %s", N_TOP, paste(top_genes, collapse = ", ")))

# ==============================================================================
# Load annotated VCC Seurat object
# ==============================================================================

message("Loading annotated VCC Seurat object...")
seurat_vcc <- qread(VCC_QS)
message(sprintf("Loaded: %d cells x %d genes", ncol(seurat_vcc), nrow(seurat_vcc)))

# ==============================================================================
# Subset to NTC cells with valid probe barcode, across all 3 lanes
# ==============================================================================

cells_ntc <- colnames(seurat_vcc)[
  !is.na(seurat_vcc$target_gene) &
  seurat_vcc$target_gene == NTC_LABEL &
  !is.na(seurat_vcc$probe_barcode)
]
message(sprintf("NTC cells with valid probe barcode: %d", length(cells_ntc)))

seurat_ntc <- subset(seurat_vcc, cells = cells_ntc)
rm(seurat_vcc); gc()

message("NTC lane distribution:")
print(table(seurat_ntc$lane))
message("NTC probe barcode distribution:")
print(table(seurat_ntc$probe_barcode))

# ==============================================================================
# Extract normalized expression for top genes
# ==============================================================================

expr <- GetAssayData(seurat_ntc, layer = "data")   # Seurat 5: layer = "data"

genes_present <- top_genes[top_genes %in% rownames(expr)]
if (length(genes_present) == 0) stop("None of the top eta² genes found in the Seurat object.")
if (length(genes_present) < length(top_genes)) {
  message("Warning: not all top genes found: ",
          paste(setdiff(top_genes, genes_present), collapse = ", "))
}

expr_sub <- as.matrix(expr[genes_present, , drop = FALSE])

# ==============================================================================
# Build long-format data frame and save
# ==============================================================================

meta <- data.frame(
  cell          = colnames(seurat_ntc),
  probe_barcode = seurat_ntc$probe_barcode,
  lane          = seurat_ntc$lane,
  stringsAsFactors = FALSE
)

df_long <- as.data.frame(t(expr_sub)) |>
  rownames_to_column("cell") |>
  pivot_longer(-cell, names_to = "gene", values_to = "expression") |>
  left_join(meta, by = "cell")

message(sprintf("Output rows: %d  (%d genes x %d cells)",
                nrow(df_long), length(genes_present), ncol(seurat_ntc)))

write.csv(df_long, OUT_CSV, row.names = FALSE)
message("Saved: ", OUT_CSV)
message("Done.")
print(sessionInfo())
