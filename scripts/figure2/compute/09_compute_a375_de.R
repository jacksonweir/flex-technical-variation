#!/usr/bin/env Rscript
# ==============================================================================
# 09_compute_a375_de.R
# FindAllMarkers (one-vs-all) across probe barcodes — A375 FLEXv2 (Figure 2B)
# ==============================================================================
#
# Each probe barcode in this dataset corresponds to a ligand treatment condition
# (PBS, ANXA1, GDF15, IFNG, MDK, POSTN, SPP1, TNFA). FindAllMarkers compares
# each condition against all others pooled, producing DE results that reflect
# both biological (ligand) and technical (barcode) effects.
#
# Parameters: Wilcoxon rank-sum, min.pct = 0.1, logfc.threshold = 0,
#             only.pos = FALSE. Significance filtering is applied at plotting
#             time (|log2FC| >= 0.3, p_val_adj < 0.05).
#
# INPUT
#   results/seurat_a375_flexv2_filtered.qs  — from 08_setup_a375_seurat.R
#
# OUTPUT
#   results/DE_A375_flexv2_FindAllMarkers_wilcox.csv
#
# USAGE
#   Rscript scripts/figure2/compute/09_compute_a375_de.R
#
# PACKAGES
#   Seurat, qs
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
})

RESULTS_DIR <- "results"
IN_QS       <- file.path(RESULTS_DIR, "seurat_a375_flexv2_filtered.qs")
OUT_CSV     <- file.path(RESULTS_DIR, "DE_A375_flexv2_FindAllMarkers_wilcox_minpct0.1.csv")

if (!file.exists(IN_QS)) {
  stop("Filtered A375 Seurat not found: ", IN_QS,
       "\nRun 08_setup_a375_seurat.R first.")
}

if (file.exists(OUT_CSV)) {
  message("Already exists, skipping: ", OUT_CSV)
  quit(status = 0)
}

# ==============================================================================
# Load and run FindAllMarkers
# ==============================================================================

message("Loading filtered A375 FLEXv2 Seurat object...")
obj <- qread(IN_QS)
message(sprintf("Cells: %d", ncol(obj)))

Idents(obj) <- obj$probe_barcode
message(sprintf("Running FindAllMarkers across %d probe barcodes (ligand conditions)...",
                length(levels(Idents(obj)))))
message("Ident levels: ", paste(levels(Idents(obj)), collapse = ", "))

df <- FindAllMarkers(obj,
  only.pos        = FALSE,
  min.pct         = 0.1,
  logfc.threshold = 0,
  test.use        = "wilcox",
  verbose         = TRUE
)
rm(obj); gc()

message(sprintf("FindAllMarkers returned %d rows", nrow(df)))

# ==============================================================================
# Save
# ==============================================================================

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
message("Saved: ", OUT_CSV)

message("\nSignificant genes (|log2FC| >= 0.3, p_adj < 0.05) per probe barcode:")
sig <- df[!is.na(df$p_val_adj) & df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.3, ]
print(sort(table(sig$cluster)))

message("Done.")
print(sessionInfo())
