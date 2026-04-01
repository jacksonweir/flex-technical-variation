# Analysis Code Reference

All code below is adapted from the paper's scripts in `scripts/`. References to the
source script are given for each section.

---

## Pairwise DE (FindMarkers)

Adapted from `scripts/figure1/compute/02_compute_pairwise_de.R`.
The actual script subsets cells by both probe barcode AND lane using `WhichCells()`,
then passes cell vectors directly to `ident.1`/`ident.2` with `group.by = NULL`.

```r
library(Seurat)

# Subset to cells of interest using WhichCells
cells_bc002 <- WhichCells(obj, expression = probe_barcode == "BC002" & lane == "1")
cells_bc011 <- WhichCells(obj, expression = probe_barcode == "BC011" & lane == "1")

de_pair <- FindMarkers(
  obj,
  ident.1         = cells_bc002,
  ident.2         = cells_bc011,
  group.by        = NULL,
  test.use        = "wilcox",
  min.pct         = 0.1,
  logfc.threshold = 0,
  verbose         = FALSE
)
de_pair$gene <- rownames(de_pair)
de_pair$p_val_adj[is.na(de_pair$p_val_adj)] <- 1

# Count significant DE genes
sum(de_pair$p_val_adj < 0.05 & abs(de_pair$avg_log2FC) >= 0.3, na.rm = TRUE)
```

Paper benchmark: BC002 vs BC011 (worst-case pair in VCC Lane 1) = **217 DE genes** at |log2FC| > 0.3, using all Lane 1 cells with Seurat 4.3.0.1.

---

## FindAllMarkers (One-vs-All Per Barcode)

Adapted from `scripts/figure1/compute/03_compute_findallmarkers_lane1.R`.

```r
library(Seurat)

# Subset to Lane 1 (or relevant lane/condition subset)
seurat_l1 <- subset(obj, subset = lane == "1" & !is.na(probe_barcode))

Idents(seurat_l1) <- seurat_l1$probe_barcode

de_all <- FindAllMarkers(
  seurat_l1,
  only.pos        = FALSE,
  min.pct         = 0.1,
  logfc.threshold = 0,
  test.use        = "wilcox",
  verbose         = TRUE
)

# Filter significant DE genes per barcode
sig <- de_all[!is.na(de_all$p_val_adj) &
              de_all$p_val_adj < 0.05 &
              abs(de_all$avg_log2FC) >= 0.3, ]
sort(table(sig$cluster))
```

---

## Eta-Squared Computation

Adapted from `scripts/figure1/compute/04_compute_variance_explained.R` (Flex v1 VCC) and
`scripts/figure2/compute/10_compute_variance_explained_a375.R` (Flex v2 A375).

Uses sparse matrix algebra — no densification required. NTC cells only for VCC.

```r
library(Seurat)
library(Matrix)
library(dplyr)

# Subset to NTC cells (VCC) or all cells (A375 Flex v2)
cells_ntc <- colnames(obj)[
  !is.na(obj$target_gene) &
  obj$target_gene == "non-targeting" &
  !is.na(obj$probe_barcode)
]
seurat_ntc <- subset(obj, cells = cells_ntc)

# Extract log-normalised expression matrix
expr <- GetAssayData(seurat_ntc, layer = "data")  # genes x cells (use slot= for Seurat v4)

# Filter to expressed genes (>= 5% of cells)
pct_expr  <- Matrix::rowMeans(expr > 0)
keep_genes <- pct_expr >= 0.05
expr_filt  <- expr[keep_genes, ]

# Eta-squared function (one-way ANOVA, sparse matrix algebra)
compute_eta_sq <- function(mat, group_vec) {
  N      <- ncol(mat)
  groups <- sort(unique(group_vec))
  K      <- length(groups)
  df1    <- K - 1
  df2    <- N - K

  grand_mean <- Matrix::rowMeans(mat)
  SS_total   <- as.numeric(Matrix::rowSums(mat^2) - N * grand_mean^2)
  SS_between <- rep(0.0, nrow(mat))

  for (g in groups) {
    idx <- which(group_vec == g)
    if (length(idx) < 2) next
    gm          <- as.numeric(Matrix::rowMeans(mat[, idx, drop = FALSE]))
    SS_between  <- SS_between + length(idx) * gm^2
  }
  SS_between <- SS_between - N * as.numeric(grand_mean)^2

  eta_sq           <- SS_between / SS_total
  eta_sq[SS_total < 1e-10] <- 0
  eta_sq[eta_sq < 0]       <- 0
  eta_sq[eta_sq > 1]       <- 1

  SS_within <- pmax(SS_total - SS_between, 0)
  F_stat    <- (SS_between / df1) / (SS_within / df2)
  F_stat[SS_total < 1e-10] <- 0

  p_val <- pf(F_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
  p_adj <- p.adjust(p_val, method = "BH")

  data.frame(gene = rownames(mat), eta_sq = eta_sq,
             F_stat = F_stat, p_val = p_val, p_adj = p_adj,
             stringsAsFactors = FALSE)
}

# Compute for probe barcode and lane
res_bc   <- compute_eta_sq(expr_filt, seurat_ntc$probe_barcode)
res_lane <- compute_eta_sq(expr_filt, seurat_ntc$lane)

cat(sprintf("Genes with eta^2_bc   >= 1%%: %d\n", sum(res_bc$eta_sq   >= 0.01)))
cat(sprintf("Genes with eta^2_lane >= 1%%: %d\n", sum(res_lane$eta_sq >= 0.01)))
```

Paper benchmarks (VCC NTC cells, all 3 lanes, all 16 probe barcodes):
- **519 genes** with η²_bc ≥ 1%
- **27 genes** with η²_lane ≥ 1%

---

## FDP Analysis

Adapted from `scripts/figure1/compute/05_compute_fdp_sweep.R`.

**Requires orthogonal barcoding** (CRISPR guide assignments in the VCC dataset).
The paper uses `TRUTH_BC = "BC005"` and `CONFOUNDED_BC = "BC008"` as the average-effect pair.

The approach:
- **Ground truth DE**: NTC cells (TRUTH_BC) vs CRISPR target cells (TRUTH_BC) — same barcode, no confounding → true positive gene set
- **Confounded DE**: CRISPR target cells (CONFOUNDED_BC) vs NTC cells (TRUTH_BC) — different barcodes, simulating a confounded design

```r
library(Seurat)

TRUTH_BC      <- "BC005"
CONFOUNDED_BC <- "BC008"
NTC_LABEL     <- "non-targeting"
LOG2FC_CUTOFFS <- seq(0, 2, by = 0.1)
PVAL_CUTOFF   <- 0.05
MIN_PCT       <- 0.1

cells_ntc <- WhichCells(obj,
  expression = target_gene == NTC_LABEL & probe_barcode == TRUTH_BC)

# For a given target gene:
tg <- "your_target_gene"

cells_truth <- WhichCells(obj,
  expression = target_gene == tg & probe_barcode == TRUTH_BC)
cells_conf  <- WhichCells(obj,
  expression = target_gene == tg & probe_barcode == CONFOUNDED_BC)

# Ground truth DE (within TRUTH_BC — no barcode confounding)
de_truth <- FindMarkers(obj, ident.1 = cells_ntc, ident.2 = cells_truth,
                        group.by = NULL, logfc.threshold = 0,
                        min.pct = MIN_PCT, test.use = "wilcox")
de_truth$gene <- rownames(de_truth)
de_truth$p_val_adj[is.na(de_truth$p_val_adj)] <- 1

# Confounded DE (target in CONFOUNDED_BC vs NTC in TRUTH_BC)
de_conf <- FindMarkers(obj, ident.1 = cells_conf, ident.2 = cells_ntc,
                       group.by = NULL, logfc.threshold = 0,
                       min.pct = MIN_PCT, test.use = "wilcox")
de_conf$gene <- rownames(de_conf)
de_conf$p_val_adj[is.na(de_conf$p_val_adj)] <- 1

# Merge and compute FDP at each log2FC threshold
df_merged <- merge(
  de_truth[, c("gene", "avg_log2FC", "p_val_adj")],
  de_conf[,  c("gene", "avg_log2FC", "p_val_adj")],
  by = "gene", suffixes = c("_truth", "_confounded")
)

calc_fdp <- function(lfc_truth, padj_truth, lfc_pred, padj_pred, lfc_cut, pval_cut) {
  gt   <- abs(lfc_truth) >= lfc_cut & padj_truth < pval_cut
  pred <- abs(lfc_pred)  >= lfc_cut & padj_pred  < pval_cut
  TP   <- sum(gt &  pred, na.rm = TRUE)
  FP   <- sum(!gt & pred, na.rm = TRUE)
  if ((TP + FP) == 0) return(NA_real_)
  FP / (TP + FP)
}

fdp_sweep <- sapply(LOG2FC_CUTOFFS, function(lfc_cut) {
  calc_fdp(df_merged$avg_log2FC_truth,      df_merged$p_val_adj_truth,
           df_merged$avg_log2FC_confounded, df_merged$p_val_adj_confounded,
           lfc_cut, PVAL_CUTOFF)
})

fdp_df <- data.frame(log2fc_cutoff = LOG2FC_CUTOFFS, fdp = fdp_sweep)
```

Script 05 loops over all target genes and saves per-gene metrics. The combined output
`results/fdp_sweep_BC005_vs_BC008/metrics_all_genes.csv` is used to compute the median FDP
shown in Figure 1H.

---

## Visualization

### Volcano Plot

Adapted from `scripts/figure1/plot/fig1bc_volcanos.R`.

```r
library(ggplot2)
library(ggthemes)
library(ggrepel)

FC_THRESH   <- 0.3
PVAL_CUTOFF <- 0.05

de_pair$log10_pval <- -log10(de_pair$p_val_adj + 1e-300)
de_pair$reg_dir    <- "Not Significant"
de_pair$reg_dir[de_pair$p_val_adj < PVAL_CUTOFF & de_pair$avg_log2FC >  FC_THRESH] <- "Upregulated"
de_pair$reg_dir[de_pair$p_val_adj < PVAL_CUTOFF & de_pair$avg_log2FC < -FC_THRESH] <- "Downregulated"

# Label top genes by combined score
sig <- de_pair$reg_dir != "Not Significant"
score     <- de_pair$log10_pval * abs(de_pair$avg_log2FC)
top_genes <- de_pair$gene[sig][order(score[sig], decreasing = TRUE)[seq_len(min(60, sum(sig)))]]
de_pair$label <- ifelse(de_pair$gene %in% top_genes, de_pair$gene, NA_character_)

ggplot(de_pair, aes(x = avg_log2FC, y = log10_pval)) +
  geom_point(aes(color = reg_dir), alpha = 0.8, size = 2, shape = 16) +
  scale_color_manual(values = c(
    "Upregulated"     = "tomato3",
    "Downregulated"   = "steelblue3",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(PVAL_CUTOFF),
             linetype = "dashed", color = "black") +
  geom_text_repel(
    aes(label = label),
    size = 6, force = 2.5, max.overlaps = 5,
    box.padding = 0.3, point.padding = 0.2, segment.size = 0.2,
    color = "black", segment.color = "black", na.rm = TRUE
  ) +
  labs(x = "Average log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_few() +
  theme(
    axis.title        = element_text(size = 20, color = "black"),
    axis.text         = element_text(size = 16, color = "black"),
    legend.position   = "none",
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )
```

### Violin Plots of Top η² Genes

Adapted from `scripts/figure1/plot/fig1d_violin_top_genes.R`.
The actual script reads from a pre-extracted CSV (`00_prep_violin_data.R`) and uses
custom ggplot2 violins split by lane. Quick diagnostic using Seurat's `VlnPlot`:

```r
VlnPlot(obj, features = c("SRP9", "SPCS1", "ISCU"),
        group.by = "probe_barcode", ncol = 3, pt.size = 0)
```

Subtle, systematic shifts across all barcodes (no discrete on/off pattern) are
characteristic of probe barcode technical variation — see Fig. 1D in the paper.

### Stacked Barplot of DE Gene Counts Per Barcode

Adapted from `scripts/figure1/plot/fig1e_stacked_barplot.R`.

```r
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

LOG2FC_CUTOFF <- 0.3
PVAL_CUTOFF   <- 0.05

counts <- de_all |>
  filter(!is.na(avg_log2FC), !is.na(p_val_adj)) |>
  filter(abs(avg_log2FC) >= LOG2FC_CUTOFF, p_val_adj < PVAL_CUTOFF) |>
  mutate(direction = case_when(
    avg_log2FC >= LOG2FC_CUTOFF  ~ "Upregulated",
    avg_log2FC <= -LOG2FC_CUTOFF ~ "Downregulated"
  )) |>
  filter(!is.na(direction)) |>
  group_by(cluster, direction) |>
  summarise(n_genes = n(), .groups = "drop") |>
  complete(cluster, direction, fill = list(n_genes = 0)) |>
  mutate(
    cluster   = factor(cluster, levels = paste0("BC", sprintf("%03d", 1:16))),
    direction = factor(direction, levels = c("Downregulated", "Upregulated"))
  )

ggplot(counts, aes(x = cluster, y = n_genes, fill = direction)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Upregulated"   = "tomato3",
                               "Downregulated" = "steelblue3")) +
  labs(x = "Probe barcode", y = "Number of DE genes", fill = NULL) +
  theme_few() +
  theme(
    axis.title        = element_text(size = 20, color = "black"),
    axis.text         = element_text(size = 16, color = "black"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right"
  )
```

### FDP Sweep Plot

Adapted from `scripts/figure1/plot/fig1h_fdp_sweep.R`.
Grey lines = individual target genes; black line = median across all genes.

```r
library(ggplot2)
library(dplyr)

# df_all has columns: target_gene, log2fc_cutoff, FDP
# (output of 05_compute_fdp_sweep.R → metrics_all_genes.csv)
df_median <- df_all |>
  group_by(log2fc_cutoff) |>
  summarise(med_fdp = median(FDP, na.rm = TRUE), .groups = "drop")

ggplot() +
  geom_line(
    data      = df_all,
    aes(x = log2fc_cutoff, y = FDP, group = target_gene),
    linewidth = 0.35, color = "grey50", alpha = 0.10
  ) +
  geom_line(
    data      = df_median,
    aes(x = log2fc_cutoff, y = med_fdp),
    linewidth = 1.6, color = "black"
  ) +
  geom_point(
    data  = df_median,
    aes(x = log2fc_cutoff, y = med_fdp),
    size  = 2.2, color = "black", alpha = 0.95
  ) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "|log2FC| cutoff", y = "False Discovery Proportion") +
  theme_classic() +
  theme(
    legend.position   = "none",
    axis.title        = element_text(size = 20, color = "black"),
    axis.text         = element_text(size = 16, color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )
```

### Cross-Dataset Reproducibility Scatter

Adapted from `scripts/figure1/plot/fig1i_cross_dataset_scatter.R`.
Compares BC002 vs BC011 log2FC between two independent datasets to show the artifact is
reproducible. Paper uses VCC (Lane 1) vs 10x PBMC 320k.

```r
library(ggplot2)
library(ggrepel)
library(dplyr)

PVAL_CUTOFF <- 0.05
FC_LINE     <- 0.3

# Merge two pairwise DE results by gene
df <- merge(
  df_dataset1[, c("gene", "avg_log2FC", "p_val_adj")],
  df_dataset2[, c("gene", "avg_log2FC", "p_val_adj")],
  by = "gene", suffixes = c("_ds1", "_ds2")
)
df$p_val_adj_ds1[is.na(df$p_val_adj_ds1)] <- 1
df$p_val_adj_ds2[is.na(df$p_val_adj_ds2)] <- 1

df$category <- with(df,
  ifelse(p_val_adj_ds1 < PVAL_CUTOFF & p_val_adj_ds2 < PVAL_CUTOFF, "Significant in both",
  ifelse(p_val_adj_ds1 < PVAL_CUTOFF, "Significant in dataset 1 only",
  ifelse(p_val_adj_ds2 < PVAL_CUTOFF, "Significant in dataset 2 only",
         "Not significant in either"))))

df$label <- ifelse(df$category != "Not significant in either", df$gene, NA_character_)

ggplot(df, aes(x = avg_log2FC_ds1, y = avg_log2FC_ds2, color = category)) +
  geom_point(alpha = 0.75, size = 1.6, shape = 16) +
  geom_vline(xintercept = c(-FC_LINE, FC_LINE),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-FC_LINE, FC_LINE),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c(
    "Significant in both"           = "red3",
    "Significant in dataset 1 only" = "deepskyblue4",
    "Significant in dataset 2 only" = "orange3",
    "Not significant in either"     = "grey"
  )) +
  geom_text_repel(
    aes(label = label),
    size = 6, force = 0.1, max.overlaps = 8,
    box.padding = 0.2, point.padding = 0.15, segment.size = 0.2,
    na.rm = TRUE, color = "black", segment.color = "black"
  ) +
  labs(x = "Average log2FC (dataset 1)", y = "Average log2FC (dataset 2)") +
  theme_classic() +
  theme(
    axis.title        = element_text(size = 20, color = "black"),
    axis.text         = element_text(size = 16, color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "none"
  )

# Spearman correlation
cor_result <- cor.test(df$avg_log2FC_ds1, df$avg_log2FC_ds2, method = "spearman")
cat(sprintf("Spearman r = %.3f (p = %.2e)\n", cor_result$estimate, cor_result$p.value))
```

Paper benchmark: VCC vs PBMC BC002 vs BC011 Spearman ρ = **0.73**, p < 0.01.

---

## Probe-Level h5 Analysis (Advanced — Python)

For investigating probe-level mechanisms using the CellRanger `molecule_info.h5` output.
Adapted from `scripts/figure2/compute/11_compute_vcc_probe_expression.R` (R) and
the Python preprocessing approach used in the analysis notebooks.

```python
import h5py
import numpy as np
from scipy.sparse import csr_matrix

CHUNK_SIZE = 50_000_000  # avoid memory overload — full file is ~155GB uncompressed

h5_path = "path/to/sample_molecule_info.h5"

with h5py.File(h5_path, "r") as f:
    # IMPORTANT: use probes/probe_id (not features/id)
    # probes/probe_id = full probe IDs like "ENSG00000000003-TSPAN6-8eab823"
    # features/id    = Ensembl gene IDs only (fewer rows, no probe-level info)
    probe_ids   = f["probes/probe_id"][:].astype(str)
    barcodes_h5 = f["barcodes"][:].astype(str)

    n_total = f["count"].shape[0]
    probe_idx_all, barcode_idx_all, count_all = [], [], []

    for start in range(0, n_total, CHUNK_SIZE):
        end = min(start + CHUNK_SIZE, n_total)
        probe_idx_all.append(f["probe_idx"][start:end])
        barcode_idx_all.append(f["barcode_idx"][start:end])
        count_all.append(f["count"][start:end])

probe_idx   = np.concatenate(probe_idx_all)
barcode_idx = np.concatenate(barcode_idx_all)
counts      = np.concatenate(count_all)

# Build sparse matrix: probes x barcodes
mat = csr_matrix((counts, (probe_idx, barcode_idx)),
                 shape=(len(probe_ids), len(barcodes_h5)))
```

**Parse probe IDs** (format: `ENSG_ID-GENE_SYMBOL-PROBE_HASH`, separator is `-` not `|`):
```python
def parse_gene_symbol(probe_id):
    parts = probe_id.split("-")
    return parts[1] if len(parts) >= 2 else probe_id

gene_symbols = [parse_gene_symbol(p) for p in probe_ids]
```

**Key finding from paper:** For SRP9, 2 probes exist — one essentially silent (0.4% detection,
0.004 mean UMIs) and one fully active (~22 UMIs/cell). The entire barcode DE effect on SRP9
is driven by the single active probe showing reduced expression in BC002. 561/1,366 target
probes (41%) have η² ≥ 1%; 469/507 genes (93%) have at least one affected probe.
