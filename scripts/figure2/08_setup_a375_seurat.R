#!/usr/bin/env Rscript
# ==============================================================================
# 08_setup_a375_seurat.R
# Filter A375 FLEXv2 Seurat object and save for downstream Figure 2 analyses
# ==============================================================================
#
# The A375 Flex v2 ligand dictionary dataset contains A375 cells treated with
# eight ligands (ANXA1, GDF15, IFNG, MDK, PBS, POSTN, SPP1, TNFA). Each ligand
# condition was assigned to a distinct Flex v2 sample barcode, which has been
# demultiplexed into a condition label stored in the `probe_barcode` metadata
# column. PBS (phosphate-buffered saline) serves as the no-treatment control.
#
# NOTE: This dataset is not publicly available. To obtain access, contact the
# Chen lab at the Broad Institute. Users must set A375_QS_PATH below to the
# location of the raw Seurat object on their system.
#
# QC filters applied:
#   nCount_RNA   > 1000
#   nFeature_RNA > 1000
#   percent.mito < 5
#
# INPUT
#   <A375_QS_PATH>  — raw A375 FLEXv2 Seurat object (user must set path)
#
# OUTPUT
#   results/seurat_a375_flexv2_filtered.qs
#
# USAGE
#   Rscript scripts/figure2/08_setup_a375_seurat.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
})

# User must set this path to the A375 FLEXv2 Seurat object on their system
A375_QS_PATH <- "/path/to/seurat_A375_flexv2.qs"

RESULTS_DIR <- "results"
OUT_QS      <- file.path(RESULTS_DIR, "seurat_a375_flexv2_filtered.qs")

MIN_COUNT   <- 1000
MIN_FEATURE <- 1000
MAX_MITO    <- 5

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(A375_QS_PATH)) {
  stop("A375 Seurat object not found: ", A375_QS_PATH,
       "\nThis dataset is not publicly available.",
       "\nContact the Chen lab (Broad Institute) for access and update A375_QS_PATH.")
}

if (file.exists(OUT_QS)) {
  message("Already exists, skipping: ", OUT_QS)
  quit(status = 0)
}

# ==============================================================================
# Load and filter
# ==============================================================================

message("Loading A375 FLEXv2 Seurat object...")
obj <- qread(A375_QS_PATH)
message(sprintf("Before filtering: %d cells", ncol(obj)))
message("Probe barcode (ligand condition) distribution:")
print(sort(table(obj$probe_barcode)))

obj_filt <- subset(obj,
  subset = nCount_RNA > MIN_COUNT &
           nFeature_RNA > MIN_FEATURE &
           percent.mito < MAX_MITO)
message(sprintf("After filtering:  %d cells (removed %d)",
                ncol(obj_filt), ncol(obj) - ncol(obj_filt)))
rm(obj); gc()

# ==============================================================================
# Save
# ==============================================================================

message("Saving filtered Seurat object...")
qsave(obj_filt, OUT_QS)
message("Saved: ", OUT_QS)

message("\nFiltered cell counts per probe barcode:")
print(sort(table(obj_filt$probe_barcode)))

message("Done.")
print(sessionInfo())
