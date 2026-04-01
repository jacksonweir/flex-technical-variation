#!/usr/bin/env Rscript
# ==============================================================================
# 11_compute_vcc_probe_expression.R
# Build probe-level expression cache for VCC NTC cells (Figures 2E)
# ==============================================================================
#
# Extracts normalised and raw count expression matrices at probe resolution
# (individual probes as features, not collapsed to gene level) from the VCC
# Lane 1 NTC dataset. These matrices are saved as a compact cache used by
# fig2e_vcc_probe_expression.R to plot per-probe expression distributions
# across probe barcodes.
#
# Probe-level data comes from the 10x cellranger output in probe-level mode,
# where each row corresponds to a unique probe ID (e.g.
# ENSG00000143742|SRP9|ac81972) rather than a gene. Only NTC cells from Lane 1
# are included, matching the VCC analysis in Figure 1.
#
# NOTE: The pre-built expression cache vcc_probe_expr_sparse_cache_10k_feat_filt.rds
# is publicly available from the Zenodo deposit (https://doi.org/10.5281/zenodo.19363777).
# Download that file and set VCC_PROBE_CACHE below to its local path. The probe-level
# Seurat object (VCC_PROBE_RDS) is an alternative source on the Chen lab servers;
# the Zenodo cache is preferred for reproducibility.
#
# IMPORTANT: Requires Matrix >= 1.7 (R 4.5+). Older Matrix versions may drop
# dgCMatrix rownames when loading the probe expression cache.
#
# INPUT  (user must set paths)
#   VCC_PROBE_CACHE — sparse expression cache from Zenodo:
#                     vcc_probe_expr_sparse_cache_10k_feat_filt.rds
#                     (rds list with $data, $counts, $barcodes fields)
#   or
#   VCC_PROBE_RDS   — probe-level Seurat object (vcc_probe_obj_lane1_ntc_10k_feat_filt.rds,
#                     Chen lab servers only)
#
# OUTPUT
#   results/vcc_probe_expression_cache.rds
#     A named list: $data (normalised expr, probes x cells),
#                   $counts (raw counts, probes x cells),
#                   $barcodes (probe_barcode per cell)
#
# USAGE
#   Rscript scripts/figure2/compute/11_compute_vcc_probe_expression.R
#
# PACKAGES
#   Seurat, Matrix
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Download from Zenodo (https://doi.org/10.5281/zenodo.19363777) and set path:
VCC_PROBE_CACHE <- "/path/to/vcc_probe_expr_sparse_cache_10k_feat_filt.rds"
# Alternatively, if you have the Chen lab probe-level Seurat object:
VCC_PROBE_RDS   <- "/path/to/vcc_probe_obj_lane1_ntc_10k_feat_filt.rds"

RESULTS_DIR <- "results"
OUT_RDS     <- file.path(RESULTS_DIR, "vcc_probe_expression_cache.rds")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

if (file.exists(OUT_RDS)) {
  message("Already exists, skipping: ", OUT_RDS)
  quit(status = 0)
}

# ==============================================================================
# Load expression data: prefer cache, fall back to Seurat object
# ==============================================================================

if (file.exists(VCC_PROBE_CACHE)) {
  message("Loading pre-built VCC probe expression cache...")
  cache    <- readRDS(VCC_PROBE_CACHE)
  expr_data   <- cache$data
  expr_counts <- cache$counts
  barcodes    <- cache$barcodes
  rm(cache); gc()
  message(sprintf("Loaded: %d probes x %d cells", nrow(expr_data), ncol(expr_data)))

} else if (file.exists(VCC_PROBE_RDS)) {
  message("Cache not found. Loading probe-level Seurat object: ", VCC_PROBE_RDS)
  obj <- readRDS(VCC_PROBE_RDS)
  message(sprintf("Object: %d probes x %d cells", nrow(obj), ncol(obj)))

  stopifnot("probe_barcode" %in% colnames(obj@meta.data))
  barcodes    <- obj$probe_barcode
  expr_data   <- as(LayerData(obj, layer = "data"),   "dgCMatrix")
  expr_counts <- as(LayerData(obj, layer = "counts"), "dgCMatrix")
  rm(obj); gc()

} else {
  stop(
    "Neither the probe-level Seurat object nor its cache was found.\n",
    "VCC_PROBE_RDS:   ", VCC_PROBE_RDS, "\n",
    "VCC_PROBE_CACHE: ", VCC_PROBE_CACHE, "\n",
    "Update these paths at the top of this script.\n",
    "Download vcc_probe_expr_sparse_cache_10k_feat_filt.rds from Zenodo:\n",
    "  https://doi.org/10.5281/zenodo.19363777\n",
    "Then set VCC_PROBE_CACHE at the top of this script."
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
