#!/usr/bin/env Rscript
# ==============================================================================
# 07_compute_pbmc_de.R
# BC002 vs BC011 DE in the PBMC 16-plex Flex v1 dataset (Figure 1I)
# ==============================================================================
#
# Runs FindMarkers() comparing probe barcode BC002 vs BC011 across all cells
# in the filtered PBMC dataset (Wilcoxon rank-sum test, min.pct=0.1,
# logfc.threshold=0). All cells (no subsetting by perturbation type, as PBMCs
# carry no guide RNA).
#
# The resulting log2FC values are compared to the equivalent VCC lane 1 result
# (from 02_compute_pairwise_de.R) to assess cross-dataset reproducibility of
# the probe barcode effect (Figure 1I).
#
# INPUT
#   results/seurat_pbmc_filtered.qs  — from 06_setup_pbmc_seurat.R
#
# OUTPUT
#   results/DE_BC002_vs_BC011_PBMC_wilcox.csv
#
# USAGE
#   Rscript scripts/figure1/compute/07_compute_pbmc_de.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
})

RESULTS_DIR <- "results"
PBMC_QS     <- file.path(RESULTS_DIR, "seurat_pbmc_filtered.qs")
OUT_CSV     <- file.path(RESULTS_DIR, "DE_BC002_vs_BC011_PBMC_wilcox.csv")
MIN_PCT     <- 0.1

if (!file.exists(PBMC_QS)) {
  stop("PBMC Seurat not found: ", PBMC_QS, "\nRun 06_setup_pbmc_seurat.R first.")
}

if (file.exists(OUT_CSV)) {
  message("Already exists, skipping: ", OUT_CSV)
  quit(status = 0)
}

# ==============================================================================
# Load
# ==============================================================================

message("Loading filtered PBMC Seurat object...")
seurat_pbmc <- qread(PBMC_QS)
message(sprintf("Loaded: %d cells", ncol(seurat_pbmc)))

# ==============================================================================
# Run DE: BC002 vs BC011
# ==============================================================================

cells_bc002 <- WhichCells(seurat_pbmc, expression = probe_barcode == "BC002")
cells_bc011 <- WhichCells(seurat_pbmc, expression = probe_barcode == "BC011")
message(sprintf("BC002: %d cells  |  BC011: %d cells", length(cells_bc002), length(cells_bc011)))

message("Running FindMarkers (Wilcoxon, min.pct=0.1)...")
de <- FindMarkers(
  seurat_pbmc,
  ident.1         = cells_bc002,
  ident.2         = cells_bc011,
  group.by        = NULL,
  test.use        = "wilcox",
  min.pct         = MIN_PCT,
  logfc.threshold = 0,
  verbose         = TRUE
)

de$gene <- rownames(de)
de$p_val_adj[is.na(de$p_val_adj)] <- 1
message(sprintf("Genes tested: %d", nrow(de)))

# ==============================================================================
# Save
# ==============================================================================

write.csv(de, OUT_CSV, row.names = FALSE)
message("Saved: ", OUT_CSV)

message("Done.")
print(sessionInfo())
