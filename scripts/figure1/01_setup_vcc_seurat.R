#!/usr/bin/env Rscript
# ==============================================================================
# 01_setup_vcc_seurat.R
# Annotate the Virtual Cell Challenge (VCC) Seurat object
# ==============================================================================
#
# Loads the merged VCC Seurat object published by the Arc Institute for the
# Virtual Cell Challenge and adds the metadata columns needed for all
# downstream analyses:
#
#   lane                    — sequencing lane (1, 2, or 3), parsed from 'batch'
#   probe_barcode_seq       — 8 bp probe barcode, extracted from cell barcode
#   probe_barcode           — probe barcode ID (BC001–BC016)
#   percent.mito            — % mitochondrial UMIs per cell
#   target_gene_probe_barcode — "<target_gene>_<probe_barcode>" convenience column
#
# INPUT
#   VCC merged Seurat object (.qs)
#   Set VCC_QS_PATH below to point to your local copy.
#
# OUTPUT
#   <RESULTS_DIR>/seurat_vcc_annotated.qs   — annotated Seurat object
#
# DATA AVAILABILITY
#   The VCC merged Seurat object is publicly available from the Arc Institute
#   Virtual Cell Challenge:
#   https://www.kaggle.com/competitions/open-problems-single-cell-perturbations
#
# USAGE
#   Rscript scripts/figure1/01_setup_vcc_seurat.R
#
# EXPECTED RUNTIME
#   ~10 minutes (dominated by reading the large .qs file)
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(tibble)
})

# ==============================================================================
# CONFIGURATION — set paths here
# ==============================================================================

# Path to the merged VCC Seurat object (.qs) downloaded from the Arc Institute
VCC_QS_PATH <- "/path/to/merged_seurat.qs"

# Directory where the annotated object will be saved
RESULTS_DIR <- "results"

# ==============================================================================
# 1. Load data
# ==============================================================================

if (!file.exists(VCC_QS_PATH)) {
  stop("VCC Seurat object not found at: ", VCC_QS_PATH,
       "\nSet VCC_QS_PATH at the top of this script to your local copy.")
}

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

message("Loading merged VCC Seurat object...")
seurat_vcc <- qread(VCC_QS_PATH)
message(sprintf("Loaded: %d cells, %d features", ncol(seurat_vcc), nrow(seurat_vcc)))

# ==============================================================================
# 2. Extract lane from batch metadata
# ==============================================================================

# Batch column format: "Flex_<lane>_<index>" e.g. "Flex_1_01" → lane = "1"
message("Extracting lane from batch metadata...")
seurat_vcc@meta.data$lane <- sub("^Flex_(\\d+)_.*$", "\\1",
                                  seurat_vcc@meta.data$batch)

message(sprintf("Lanes found: %s",
                paste(sort(unique(seurat_vcc@meta.data$lane)), collapse = ", ")))

# ==============================================================================
# 3. Extract probe barcode sequence from cell barcode
# ==============================================================================

# Cell barcodes have the format: <16 bp barcode>-<suffix>
# The probe barcode sequence is the last 8 bp of the 16 bp barcode portion.
message("Extracting probe barcode sequences from cell barcodes...")

bc_16bp <- sub("-.*$", "", rownames(seurat_vcc@meta.data))
seurat_vcc@meta.data$probe_barcode_seq <- substr(bc_16bp,
                                                  nchar(bc_16bp) - 7,
                                                  nchar(bc_16bp))

# ==============================================================================
# 4. Map 8 bp sequence to probe barcode ID (BC001–BC016)
# ==============================================================================

# Reference table: 8 bp sequences for each of the 16 Flex v1 probe barcodes
# Source: 10x Genomics probe-barcodes-fixed-rna-profiling.txt
probe_barcode_ref <- tibble(
  probe_barcode = paste0("BC", sprintf("%03d", 1:16)),
  sequence      = c("ACTTTAGG", "AACGGGAA", "AGTAGGCT", "ATGTTGAC",
                    "ACAGACCT", "ATCCCAAC", "AAGTAGAG", "AGCTGTGA",
                    "ACAGTCTG", "AGTGAGTG", "AGAGGCAA", "ACTACTCA",
                    "ATACGTCA", "ATCATGTG", "AACGCCGA", "ATTCGGTT")
)

message("Mapping sequences to probe barcode IDs...")
seurat_vcc@meta.data$probe_barcode <- probe_barcode_ref$probe_barcode[
  match(seurat_vcc@meta.data$probe_barcode_seq, probe_barcode_ref$sequence)
]

n_matched   <- sum(!is.na(seurat_vcc@meta.data$probe_barcode))
n_unmatched <- sum(is.na(seurat_vcc@meta.data$probe_barcode))
message(sprintf("Matched: %d cells  |  Unmatched: %d cells", n_matched, n_unmatched))

if (n_unmatched > 0) {
  warning(sprintf("%d cells could not be matched to a probe barcode ID.", n_unmatched))
}

# ==============================================================================
# 5. Calculate mitochondrial gene percentage
# ==============================================================================

message("Calculating mitochondrial gene percentages...")
mito_genes <- grep("^MT-", rownames(seurat_vcc), value = TRUE)
message(sprintf("Mitochondrial genes found: %d", length(mito_genes)))

seurat_vcc[["percent.mito"]] <- PercentageFeatureSet(seurat_vcc, features = mito_genes)

# ==============================================================================
# 6. Create combined target_gene + probe_barcode column
# ==============================================================================

seurat_vcc@meta.data$target_gene_probe_barcode <- paste0(
  seurat_vcc@meta.data$target_gene, "_", seurat_vcc@meta.data$probe_barcode
)

# ==============================================================================
# 7. Validation summary
# ==============================================================================

message("\n=== Annotation complete ===")
message(sprintf("Cells: %d  |  Features: %d", ncol(seurat_vcc), nrow(seurat_vcc)))
message(sprintf("Lanes: %s",
                paste(sort(unique(seurat_vcc@meta.data$lane)), collapse = ", ")))
message(sprintf("Probe barcodes: %s",
                paste(sort(unique(na.omit(seurat_vcc@meta.data$probe_barcode))),
                      collapse = ", ")))
message(sprintf("Cells per probe barcode:"))
print(table(seurat_vcc@meta.data$probe_barcode))

message(sprintf("\nMitochondrial %%: min=%.2f  median=%.2f  max=%.2f",
                min(seurat_vcc@meta.data$percent.mito,    na.rm = TRUE),
                median(seurat_vcc@meta.data$percent.mito, na.rm = TRUE),
                max(seurat_vcc@meta.data$percent.mito,    na.rm = TRUE)))

# ==============================================================================
# 8. Save annotated object
# ==============================================================================

out_path <- file.path(RESULTS_DIR, "seurat_vcc_annotated.qs")
message(sprintf("\nSaving annotated Seurat object to:\n  %s", out_path))
qsave(seurat_vcc, out_path)
message("Done.")

print(sessionInfo())
