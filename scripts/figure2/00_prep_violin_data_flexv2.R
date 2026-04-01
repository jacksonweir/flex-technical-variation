#!/usr/bin/env Rscript
# ==============================================================================
# 00_prep_violin_data_flexv2.R
# Extract per-cell normalized expression of SRP9 and SPCS1 from the filtered
# A375 Flex v2 Seurat object for Figure 2C violin plots.
#
# These are the same genes plotted in Figure 1D (VCC Flex v1). Plotting them
# in A375 Flex v2 allows direct comparison of the probe barcode artifact
# magnitude between kit versions.
#
# Run with mel_spatial conda env (has qs + Seurat 5):
#   conda run -n mel_spatial Rscript scripts/figure2/00_prep_violin_data_flexv2.R
#
# INPUT
#   results/seurat_a375_flexv2_filtered.qs
#     — from 08_setup_a375_seurat.R, OR place the Zenodo file
#       seurat_A375_flexv2_filtered.qs directly in results/
#       (renaming from upper A375 to lower a375 in the filename)
#
# OUTPUT
#   results/a375_flexv2_violin_data.csv
#     long-format: cell, gene, expression, probe_barcode
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
QS_PATH     <- file.path(RESULTS_DIR, "seurat_a375_flexv2_filtered.qs")
OUT_CSV     <- file.path(RESULTS_DIR, "a375_flexv2_violin_data.csv")

GENES <- c("SRP9", "SPCS1")

if (!file.exists(QS_PATH)) {
  stop("A375 Flex v2 Seurat not found: ", QS_PATH,
       "\nEither:\n",
       "  (a) Run scripts/figure2/08_setup_a375_seurat.R, or\n",
       "  (b) Download seurat_A375_flexv2_filtered.qs from Zenodo\n",
       "      (https://zenodo.org/records/19363777) and place it in results/\n",
       "      as seurat_a375_flexv2_filtered.qs")
}

# ==============================================================================
# Load Seurat object
# ==============================================================================

message("Loading filtered A375 Flex v2 Seurat object...")
obj <- qread(QS_PATH)
message(sprintf("Loaded: %d cells x %d genes", ncol(obj), nrow(obj)))

message("Probe barcode distribution:")
print(table(obj$probe_barcode))

# ==============================================================================
# Check genes present
# ==============================================================================

genes_present <- GENES[GENES %in% rownames(obj)]
missing       <- setdiff(GENES, genes_present)
if (length(missing) > 0) {
  message("Warning: genes not found in object: ", paste(missing, collapse = ", "))
}
if (length(genes_present) == 0) stop("None of the target genes found in object.")
message("Extracting expression for: ", paste(genes_present, collapse = ", "))

# ==============================================================================
# Extract normalized expression and save long-format CSV
# ==============================================================================

expr     <- GetAssayData(obj, layer = "data")   # Seurat 5: layer = "data"
expr_sub <- as.matrix(expr[genes_present, , drop = FALSE])

meta <- data.frame(
  cell          = colnames(obj),
  probe_barcode = obj$probe_barcode,
  stringsAsFactors = FALSE
)

df_long <- as.data.frame(t(expr_sub)) |>
  rownames_to_column("cell") |>
  pivot_longer(-cell, names_to = "gene", values_to = "expression") |>
  left_join(meta, by = "cell")

message(sprintf("Output rows: %d  (%d genes x %d cells)",
                nrow(df_long), length(genes_present), ncol(obj)))

write.csv(df_long, OUT_CSV, row.names = FALSE)
message("Saved: ", OUT_CSV)
message("Done.")
print(sessionInfo())
