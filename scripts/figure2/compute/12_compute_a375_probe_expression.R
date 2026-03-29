#!/usr/bin/env Rscript
# ==============================================================================
# 12_compute_a375_probe_expression.R
# Build probe-level expression cache for A375 FLEXv2 PBS cells (Figure 2F)
# ==============================================================================
#
# Extracts normalised and raw count expression matrices at probe resolution
# from the A375 FLEXv2 PBS-treated cell dataset. Only PBS cells are included —
# these are the no-treatment control cells in which all probe barcodes should
# represent biologically equivalent cells, allowing probe barcode technical
# effects to be assessed in isolation (analogous to NTC cells in the VCC
# experiment).
#
# The probe-level Seurat object (`flex_v2_probe_obj_PBS_10k_feat_filt.rds`)
# uses `probe_barcode_number` as the sample barcode identifier — a numeric
# label (e.g. "A-A01", "A-A10") corresponding to Flex v2 sample barcodes.
#
# NOTE: This data is stored on the Chen lab servers and is not publicly
# available. Set A375_PROBE_RDS and A375_PROBE_CACHE below to the correct
# paths on your system.
#
# INPUT  (user must set paths)
#   A375_PROBE_RDS   — probe-level Seurat (flex_v2_probe_obj_PBS_10k_feat_filt.rds)
#   or
#   A375_PROBE_CACHE — pre-built sparse expression cache (rds list with $data,
#                      $counts, $barcodes fields)
#
# OUTPUT
#   results/a375_probe_expression_cache.rds
#     A named list: $data (normalised expr, probes x cells),
#                   $counts (raw counts, probes x cells),
#                   $barcodes (probe_barcode_number per cell)
#
# USAGE
#   Rscript scripts/figure2/compute/12_compute_a375_probe_expression.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Set these paths to the data on your system (Chen lab servers)
A375_PROBE_RDS   <- "/path/to/flex_v2_probe_obj_PBS_10k_feat_filt.rds"
A375_PROBE_CACHE <- "/path/to/flex_v2_probe_expr_sparse_cache_10k_feat_filt.rds"

RESULTS_DIR <- "results"
OUT_RDS     <- file.path(RESULTS_DIR, "a375_probe_expression_cache.rds")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

if (file.exists(OUT_RDS)) {
  message("Already exists, skipping: ", OUT_RDS)
  quit(status = 0)
}

# ==============================================================================
# Load expression data: prefer cache, fall back to Seurat object
# ==============================================================================

if (file.exists(A375_PROBE_CACHE)) {
  message("Loading pre-built A375 probe expression cache...")
  cache       <- readRDS(A375_PROBE_CACHE)
  expr_data   <- cache$data
  expr_counts <- cache$counts
  barcodes    <- cache$barcodes
  rm(cache); gc()
  message(sprintf("Loaded: %d probes x %d cells", nrow(expr_data), ncol(expr_data)))

} else if (file.exists(A375_PROBE_RDS)) {
  message("Cache not found. Loading probe-level Seurat object: ", A375_PROBE_RDS)
  obj <- readRDS(A375_PROBE_RDS)
  message(sprintf("Object: %d probes x %d cells", nrow(obj), ncol(obj)))

  stopifnot("probe_barcode_number" %in% colnames(obj@meta.data))
  barcodes    <- obj$probe_barcode_number
  expr_data   <- as(LayerData(obj, layer = "data"),   "dgCMatrix")
  expr_counts <- as(LayerData(obj, layer = "counts"), "dgCMatrix")
  rm(obj); gc()

} else {
  stop(
    "Neither the probe-level Seurat object nor its cache was found.\n",
    "A375_PROBE_RDS:   ", A375_PROBE_RDS, "\n",
    "A375_PROBE_CACHE: ", A375_PROBE_CACHE, "\n",
    "Update these paths at the top of this script.\n",
    "This data is stored on the Chen lab servers and is not publicly available."
  )
}

# ==============================================================================
# Save cache to results/
# ==============================================================================

message("Saving probe expression cache to: ", OUT_RDS)
saveRDS(list(data = expr_data, counts = expr_counts, barcodes = barcodes),
        file = OUT_RDS)
message(sprintf("Saved: %d probes x %d cells, %d unique probe barcodes",
                nrow(expr_data), ncol(expr_data), length(unique(barcodes))))

message("Done.")
print(sessionInfo())
