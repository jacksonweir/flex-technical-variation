#!/usr/bin/env Rscript
# ==============================================================================
# 13_compute_probe_histogram_data.R
# Compute per-gene probe counts stratified by DE status for Figure 2D.
#
# For each expressed gene in VCC (Flex v1) and A375 (Flex v2), counts how many
# unique probes target it in the respective probe set manifests, and classifies
# genes as DE or not DE based on a FindAllMarkers analysis across probe barcodes.
#
# VCC: "DE" = significant (p_val_adj < 0.05, |avg_log2FC| >= 0.3) in a
#   one-vs-all FindAllMarkers across all 16 probe barcodes in NTC cells
#   (all 3 lanes combined). Expressed = detected in >= 10% of NTC cells.
#
# A375: "DE" = significant in the FindAllMarkers from 09_compute_a375_de.R
#   (all probe barcodes, min.pct = 0.1). Expressed = detected in >= 10%
#   of cells (pct_expressed >= 0.10 from 10_compute_variance_explained_a375.R).
#
# PREREQUISITES
#   results/seurat_vcc_annotated.qs           — 01_setup_vcc_seurat.R
#   results/DE_A375_flexv2_FindAllMarkers_wilcox_minpct0.1.csv — 09_compute_a375_de.R
#   results/variance_explained_probe_barcode_A375_flexv2.csv   — 10_compute_variance_explained_a375.R
#   data/vcc_probe_manifest.csv               — Flex v1 probe set (10x download; see data/README.md)
#   data/a375_v2_probe_manifest.csv           — Flex v2 probe set (10x download; see data/README.md)
#
# OUTPUT
#   results/vcc_genes_min_0.1.csv             — VCC NTC genes detected >= 10%
#   results/DE_vcc_NTC_allLanes_FindAllMarkers_wilcox.csv  — full FindAllMarkers result
#   results/probe_counts_vcc.csv              — per-gene probe counts (VCC, Flex v1)
#   results/probe_counts_a375.csv             — per-gene probe counts (A375, Flex v2)
#     columns: gene_name, n_probe_ids, is_de
#
# USAGE
#   Rscript scripts/figure2/compute/13_compute_probe_histogram_data.R
#
# PACKAGES
#   Seurat, qs, Matrix, dplyr
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(Matrix)
  library(dplyr)
})

RESULTS_DIR   <- "results"
VCC_QS        <- file.path(RESULTS_DIR, "seurat_vcc_annotated.qs")
A375_DE_CSV   <- file.path(RESULTS_DIR, "DE_A375_flexv2_FindAllMarkers_wilcox_minpct0.1.csv")
A375_ETA_CSV  <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_A375_flexv2.csv")
VCC_MANIFEST  <- "data/vcc_probe_manifest.csv"
A375_MANIFEST <- "data/a375_v2_probe_manifest.csv"

OUT_VCC_GENES  <- file.path(RESULTS_DIR, "vcc_genes_min_0.1.csv")
OUT_VCC_DE     <- file.path(RESULTS_DIR, "DE_vcc_NTC_allLanes_FindAllMarkers_wilcox.csv")
OUT_VCC_COUNTS <- file.path(RESULTS_DIR, "probe_counts_vcc.csv")
OUT_A375_COUNTS <- file.path(RESULTS_DIR, "probe_counts_a375.csv")

NTC_LABEL     <- "non-targeting"
MIN_PCT       <- 0.1
MIN_PCT_EXPR  <- 0.10
PVAL_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 0.3

for (f in c(VCC_QS, A375_DE_CSV, A375_ETA_CSV, VCC_MANIFEST, A375_MANIFEST)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f,
         "\nSee data/README.md for download instructions and script prerequisites.")
  }
}

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. VCC: load NTC cells and compute expressed genes (>= 10% detection)
# ==============================================================================

message("=== VCC (Flex v1) ===")
message("Loading annotated VCC Seurat object...")
seurat_vcc <- qread(VCC_QS)
message(sprintf("Loaded: %d cells", ncol(seurat_vcc)))

cells_ntc <- colnames(seurat_vcc)[
  !is.na(seurat_vcc$target_gene) &
  seurat_vcc$target_gene == NTC_LABEL &
  !is.na(seurat_vcc$probe_barcode)
]
message(sprintf("NTC cells (all 3 lanes): %d", length(cells_ntc)))
seurat_ntc <- subset(seurat_vcc, cells = cells_ntc)
rm(seurat_vcc); gc()

message("NTC cells per lane:")
print(table(seurat_ntc$lane))
message("NTC cells per probe barcode:")
print(table(seurat_ntc$probe_barcode))

# Expressed genes: >= 10% detection in NTC cells
expr_ntc    <- GetAssayData(seurat_ntc, layer = "counts")
pct_detect  <- Matrix::rowMeans(expr_ntc > 0)
expressed_genes <- names(pct_detect[pct_detect >= MIN_PCT_EXPR])
message(sprintf("Genes detected in >= %.0f%% of NTC cells: %d / %d",
                MIN_PCT_EXPR * 100, length(expressed_genes), nrow(seurat_ntc)))

write.csv(data.frame(gene = expressed_genes), OUT_VCC_GENES, row.names = FALSE)
message("Saved: ", OUT_VCC_GENES)

# ==============================================================================
# 2. VCC: run FindAllMarkers across all 16 probe barcodes in NTC cells
# ==============================================================================

if (file.exists(OUT_VCC_DE)) {
  message("FindAllMarkers result already exists, loading: ", OUT_VCC_DE)
  de_vcc <- read.csv(OUT_VCC_DE)
} else {
  message("Running FindAllMarkers (NTC cells, all 3 lanes, probe_barcode as identity)...")
  message("  Parameters: Wilcoxon, min.pct = ", MIN_PCT, ", logfc.threshold = 0")
  Idents(seurat_ntc) <- seurat_ntc$probe_barcode
  de_vcc <- FindAllMarkers(
    seurat_ntc,
    only.pos        = FALSE,
    min.pct         = MIN_PCT,
    logfc.threshold = 0,
    test.use        = "wilcox",
    verbose         = TRUE
  )
  write.csv(de_vcc, OUT_VCC_DE, row.names = FALSE)
  message("Saved: ", OUT_VCC_DE)
}
rm(seurat_ntc); gc()

de_genes_vcc <- unique(de_vcc$gene[
  !is.na(de_vcc$p_val_adj) &
  de_vcc$p_val_adj < PVAL_CUTOFF &
  abs(de_vcc$avg_log2FC) >= LOG2FC_CUTOFF
])
message(sprintf("VCC DE genes (|log2FC| >= %.1f, p_adj < %.2f): %d",
                LOG2FC_CUTOFF, PVAL_CUTOFF, length(de_genes_vcc)))

# ==============================================================================
# 3. VCC: probe counts per expressed gene
# ==============================================================================

message("Reading Flex v1 probe manifest...")
manifest_vcc     <- read.csv(VCC_MANIFEST, comment.char = "#")
manifest_vcc_inc <- manifest_vcc[manifest_vcc$included == TRUE, ]
manifest_vcc_inc$gene_name <- sapply(strsplit(manifest_vcc_inc$probe_id, "\\|"), `[`, 2)

probe_counts_vcc <- manifest_vcc_inc |>
  group_by(gene_name) |>
  summarise(n_probe_ids = n_distinct(probe_id), .groups = "drop") |>
  filter(gene_name %in% expressed_genes) |>
  mutate(is_de = ifelse(gene_name %in% de_genes_vcc, "DE", "Not DE"))

message(sprintf("VCC: %d expressed genes | DE: %d | Not DE: %d",
  nrow(probe_counts_vcc),
  sum(probe_counts_vcc$is_de == "DE"),
  sum(probe_counts_vcc$is_de == "Not DE")))

write.csv(probe_counts_vcc, OUT_VCC_COUNTS, row.names = FALSE)
message("Saved: ", OUT_VCC_COUNTS)

# ==============================================================================
# 4. A375 (Flex v2): load DE genes from script 09, expressed genes from script 10
# ==============================================================================

message("\n=== A375 (Flex v2) ===")

de_a375     <- read.csv(A375_DE_CSV)
de_genes_a375 <- unique(de_a375$gene[
  !is.na(de_a375$p_val_adj) &
  de_a375$p_val_adj < PVAL_CUTOFF &
  abs(de_a375$avg_log2FC) >= LOG2FC_CUTOFF
])
message(sprintf("A375 DE genes (|log2FC| >= %.1f, p_adj < %.2f): %d",
                LOG2FC_CUTOFF, PVAL_CUTOFF, length(de_genes_a375)))

eta_a375          <- read.csv(A375_ETA_CSV)
expressed_a375    <- eta_a375$gene[eta_a375$pct_expressed >= MIN_PCT_EXPR]
message(sprintf("A375 expressed genes (>= %.0f%%): %d",
                MIN_PCT_EXPR * 100, length(expressed_a375)))

# ==============================================================================
# 5. A375: probe counts per expressed gene
# ==============================================================================

message("Reading Flex v2 probe manifest...")
manifest_a375     <- read.csv(A375_MANIFEST, comment.char = "#")
manifest_a375_inc <- manifest_a375[manifest_a375$included == TRUE, ]

probe_counts_a375 <- manifest_a375_inc |>
  group_by(gene_name) |>
  summarise(n_probe_ids = n_distinct(probe_id), .groups = "drop") |>
  filter(gene_name %in% expressed_a375) |>
  mutate(is_de = ifelse(gene_name %in% de_genes_a375, "DE", "Not DE"))

message(sprintf("A375: %d expressed genes | DE: %d | Not DE: %d",
  nrow(probe_counts_a375),
  sum(probe_counts_a375$is_de == "DE"),
  sum(probe_counts_a375$is_de == "Not DE")))

write.csv(probe_counts_a375, OUT_A375_COUNTS, row.names = FALSE)
message("Saved: ", OUT_A375_COUNTS)

message("\nDone.")
print(sessionInfo())
