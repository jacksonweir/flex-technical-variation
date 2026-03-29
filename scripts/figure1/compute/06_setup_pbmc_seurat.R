#!/usr/bin/env Rscript
# ==============================================================================
# 06_setup_pbmc_seurat.R
# Load and annotate the public 10x PBMC 16-plex Flex v1 dataset
# ==============================================================================
#
# Loads the CellRanger output for the publicly available 320k Human PBMC 16-plex
# GEM-X Flex v1 dataset from 10x Genomics, creates a Seurat object, and adds
# the following metadata columns:
#
#   probe_barcode_seq   8 bp probe barcode sequence from the last 8 bases of
#                       the 16 bp cell barcode (same extraction logic as VCC)
#   probe_barcode       probe barcode ID (BC001–BC016), mapped from sequence
#   percent.mito        percent mitochondrial UMIs per cell
#
# Cells are QC-filtered: nFeature_RNA > 200 and percent.mito < 5%.
# The filtered object is log-normalised (NormalizeData defaults).
#
# DATA
#   The dataset is publicly available from 10x Genomics:
#   "320k Human PBMCs, 16-plex GEM-X Flex, Gene Expression"
#   Download the filtered_feature_bc_matrix.h5 file and set PBMC_H5_PATH below.
#
# INPUT
#   PBMC_H5_PATH  — filtered_feature_bc_matrix.h5 (set below)
#
# OUTPUT
#   results/seurat_pbmc_filtered.qs
#
# USAGE
#   Rscript scripts/figure1/compute/06_setup_pbmc_seurat.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(tibble)
})

# ==============================================================================
# CONFIGURATION — set path to downloaded PBMC h5 file
# ==============================================================================

PBMC_H5_PATH <- "/path/to/filtered_feature_bc_matrix.h5"

RESULTS_DIR  <- "results"
OUT_QS       <- file.path(RESULTS_DIR, "seurat_pbmc_filtered.qs")

MIN_FEATURES <- 200
MAX_PCT_MITO <- 5

# ==============================================================================
# Probe barcode reference (same as VCC)
# ==============================================================================

probe_barcode_ref <- tibble(
  probe_barcode = paste0("BC", sprintf("%03d", 1:16)),
  sequence      = c("ACTTTAGG", "AACGGGAA", "AGTAGGCT", "ATGTTGAC",
                    "ACAGACCT", "ATCCCAAC", "AAGTAGAG", "AGCTGTGA",
                    "ACAGTCTG", "AGTGAGTG", "AGAGGCAA", "ACTACTCA",
                    "ATACGTCA", "ATCATGTG", "AACGCCGA", "ATTCGGTT")
)

# ==============================================================================
# Load
# ==============================================================================

if (!file.exists(PBMC_H5_PATH)) {
  stop("PBMC h5 not found: ", PBMC_H5_PATH,
       "\nSet PBMC_H5_PATH at the top of this script.")
}

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

message("Loading PBMC h5 file...")
counts <- Read10X_h5(PBMC_H5_PATH)
seurat_pbmc <- CreateSeuratObject(counts = counts, project = "PBMC_16plex_Flexv1")
message(sprintf("Loaded: %d cells x %d features", ncol(seurat_pbmc), nrow(seurat_pbmc)))

# ==============================================================================
# Add probe barcode (last 8 bp of 16 bp cell barcode)
# ==============================================================================

message("Extracting probe barcode sequences...")
bc_16bp <- sub("-.*$", "", colnames(seurat_pbmc))
seurat_pbmc$probe_barcode_seq <- substr(bc_16bp, nchar(bc_16bp) - 7, nchar(bc_16bp))

seurat_pbmc$probe_barcode <- probe_barcode_ref$probe_barcode[
  match(seurat_pbmc$probe_barcode_seq, probe_barcode_ref$sequence)
]

n_matched   <- sum(!is.na(seurat_pbmc$probe_barcode))
n_unmatched <- sum(is.na(seurat_pbmc$probe_barcode))
message(sprintf("Matched: %d cells  |  Unmatched: %d cells", n_matched, n_unmatched))

if (n_unmatched > 0) {
  warning(sprintf("%d cells could not be matched to a probe barcode.", n_unmatched))
}

# ==============================================================================
# QC metrics and filtering
# ==============================================================================

message("Calculating mitochondrial gene percentages...")
mito_genes <- grep("^MT-", rownames(seurat_pbmc), value = TRUE)
seurat_pbmc[["percent.mito"]] <- PercentageFeatureSet(seurat_pbmc, features = mito_genes)

message(sprintf("Before filtering: %d cells", ncol(seurat_pbmc)))
seurat_pbmc <- subset(seurat_pbmc,
  subset = nFeature_RNA > MIN_FEATURES & percent.mito < MAX_PCT_MITO)
message(sprintf("After filtering (nFeature > %d, pct.mito < %d%%): %d cells",
                MIN_FEATURES, MAX_PCT_MITO, ncol(seurat_pbmc)))

# ==============================================================================
# Normalize
# ==============================================================================

message("Normalizing...")
seurat_pbmc <- NormalizeData(seurat_pbmc, verbose = FALSE)

# ==============================================================================
# Summary and save
# ==============================================================================

message("\nCells per probe barcode:")
print(table(seurat_pbmc$probe_barcode))

qsave(seurat_pbmc, OUT_QS)
message("Saved: ", OUT_QS)

message("Done.")
print(sessionInfo())
