#!/usr/bin/env Rscript
# ==============================================================================
# 03_compute_findallmarkers_lane1.R
# One-vs-all DE for each probe barcode in Lane 1 (Figure 1E)
# ==============================================================================
#
# Subsets the VCC dataset to Lane 1 cells (all perturbations included), sets
# probe barcode as the identity, and runs FindAllMarkers() with the Wilcoxon
# rank-sum test. Each probe barcode is compared against all other probe barcodes
# pooled together.
#
# The resulting CSV is used by plot/fig1e_stacked_barplot.R to count DE genes
# per probe barcode after applying post-hoc |log2FC| >= 0.3 and
# p_val_adj < 0.05 filters.
#
# INPUT
#   results/seurat_vcc_annotated.qs  — from 01_setup_vcc_seurat.R
#
# OUTPUT
#   results/DE_FindAllMarkers_lane1_wilcox.csv
#
# USAGE
#   Rscript scripts/figure1/compute/03_compute_findallmarkers_lane1.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
})

RESULTS_DIR <- "results"
VCC_QS      <- file.path(RESULTS_DIR, "seurat_vcc_annotated.qs")
OUT_CSV     <- file.path(RESULTS_DIR, "DE_FindAllMarkers_lane1_wilcox.csv")

if (!file.exists(VCC_QS)) {
  stop("VCC Seurat not found: ", VCC_QS, "\nRun 01_setup_vcc_seurat.R first.")
}

if (file.exists(OUT_CSV)) {
  message("Already exists, skipping: ", OUT_CSV)
  quit(status = 0)
}

# ==============================================================================
# Load and subset to Lane 1
# ==============================================================================

message("Loading annotated VCC Seurat object...")
seurat_vcc <- qread(VCC_QS)
message(sprintf("Loaded: %d cells", ncol(seurat_vcc)))

seurat_l1 <- subset(seurat_vcc, subset = lane == "1" & !is.na(probe_barcode))
rm(seurat_vcc); gc()

message(sprintf("Lane 1 subset: %d cells across %d probe barcodes",
                ncol(seurat_l1),
                length(unique(seurat_l1$probe_barcode))))
print(table(seurat_l1$probe_barcode))

# ==============================================================================
# Run FindAllMarkers (one vs all, each probe barcode)
# ==============================================================================

Idents(seurat_l1) <- seurat_l1$probe_barcode

message("Running FindAllMarkers (Wilcoxon, min.pct=0.1)...")
de_all <- FindAllMarkers(
  seurat_l1,
  only.pos        = FALSE,
  min.pct         = 0.1,
  logfc.threshold = 0,
  test.use        = "wilcox",
  verbose         = TRUE
)

message(sprintf("Total rows: %d", nrow(de_all)))

# ==============================================================================
# Save
# ==============================================================================

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(de_all, OUT_CSV, row.names = FALSE)
message("Saved: ", OUT_CSV)

message("Done.")
print(sessionInfo())
