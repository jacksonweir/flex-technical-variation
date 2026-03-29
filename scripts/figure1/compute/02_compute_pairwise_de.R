#!/usr/bin/env Rscript
# ==============================================================================
# 02_compute_pairwise_de.R
# Compute two pairwise DE comparisons from the VCC dataset
# ==============================================================================
#
# Runs two FindMarkers() calls (Wilcoxon, min.pct=0.1, logfc.threshold=0):
#
#   Comparison 1 — Lane 1 vs Lane 2 restricted to BC001 cells only
#     Isolates the lane effect within a single probe barcode.
#     → results/DE_lane1_vs_lane2_BC001_wilcox.csv    (Figure 1B)
#
#   Comparison 2 — BC002 vs BC011 restricted to Lane 1 cells only
#     Isolates the probe barcode effect within a single sequencing lane.
#     → results/DE_BC002_vs_BC011_lane1_wilcox.csv    (Figure 1C, 1I)
#
# Both outputs are cached: if the CSV already exists the comparison is skipped.
#
# INPUT
#   results/seurat_vcc_annotated.qs  — from 01_setup_vcc_seurat.R
#
# OUTPUT
#   results/DE_lane1_vs_lane2_BC001_wilcox.csv
#   results/DE_BC002_vs_BC011_lane1_wilcox.csv
#
# USAGE
#   Rscript scripts/figure1/compute/02_compute_pairwise_de.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
})

RESULTS_DIR <- "results"
VCC_QS      <- file.path(RESULTS_DIR, "seurat_vcc_annotated.qs")
MIN_PCT     <- 0.1

OUT_LANE <- file.path(RESULTS_DIR, "DE_lane1_vs_lane2_BC001_wilcox.csv")
OUT_BC   <- file.path(RESULTS_DIR, "DE_BC002_vs_BC011_lane1_wilcox.csv")

if (!file.exists(VCC_QS)) {
  stop("VCC Seurat not found: ", VCC_QS, "\nRun 01_setup_vcc_seurat.R first.")
}

# ==============================================================================
# Load
# ==============================================================================

message("Loading annotated VCC Seurat object...")
seurat_vcc <- qread(VCC_QS)
message(sprintf("Loaded: %d cells", ncol(seurat_vcc)))

run_de <- function(seurat, cells1, cells2, label) {
  message(sprintf("Running DE: %s  (n=%d vs n=%d)...", label, length(cells1), length(cells2)))
  de <- FindMarkers(
    seurat,
    ident.1         = cells1,
    ident.2         = cells2,
    group.by        = NULL,
    test.use        = "wilcox",
    min.pct         = MIN_PCT,
    logfc.threshold = 0,
    verbose         = FALSE
  )
  de$gene <- rownames(de)
  de$p_val_adj[is.na(de$p_val_adj)] <- 1
  message(sprintf("  %d genes tested", nrow(de)))
  de
}

# ==============================================================================
# 1. Lane 1 vs Lane 2 — BC001 cells only (Figure 1B)
# ==============================================================================

if (file.exists(OUT_LANE)) {
  message("Already exists, skipping: ", OUT_LANE)
} else {
  cells_l1 <- WhichCells(seurat_vcc, expression = probe_barcode == "BC001" & lane == "1")
  cells_l2 <- WhichCells(seurat_vcc, expression = probe_barcode == "BC001" & lane == "2")
  de_lane  <- run_de(seurat_vcc, cells_l1, cells_l2, "Lane1 vs Lane2 (BC001 only)")
  write.csv(de_lane, OUT_LANE, row.names = FALSE)
  message("Saved: ", OUT_LANE)
}

# ==============================================================================
# 2. BC002 vs BC011 — Lane 1 only (Figure 1C, 1I)
# ==============================================================================

if (file.exists(OUT_BC)) {
  message("Already exists, skipping: ", OUT_BC)
} else {
  cells_bc002 <- WhichCells(seurat_vcc, expression = probe_barcode == "BC002" & lane == "1")
  cells_bc011 <- WhichCells(seurat_vcc, expression = probe_barcode == "BC011" & lane == "1")
  de_bc       <- run_de(seurat_vcc, cells_bc002, cells_bc011, "BC002 vs BC011 (Lane 1 only)")
  write.csv(de_bc, OUT_BC, row.names = FALSE)
  message("Saved: ", OUT_BC)
}

message("Done.")
print(sessionInfo())
