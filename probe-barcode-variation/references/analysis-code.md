# Analysis Code Reference

## Pairwise DE (FindMarkers)

For focused comparison between a specific barcode pair:

```r
de_pair <- FindMarkers(
  obj,
  ident.1         = "BC002",
  ident.2         = "BC011",
  min.pct         = 0.1,
  logfc.threshold = 0,
  test.use        = "wilcox"
)
de_pair$gene <- rownames(de_pair)

# Count significant DE genes
sum(de_pair$p_val_adj < 0.05 & abs(de_pair$avg_log2FC) >= 0.3, na.rm = TRUE)
```

Paper benchmark: BC002 vs BC011 (worst-case pair in VCC Lane 1) = **217 DE genes** at |log2FC| > 0.3, using all cells with Seurat 4.3.0.1.

---

## Eta-Squared Computation

Eta-squared (η²) = SS_between / SS_total from one-way ANOVA across probe barcode groups.
Uses sparse matrix algebra — no matrix densification required. Matches paper scripts `04_compute_variance_explained.R` and `10_compute_variance_explained_a375.R`.

```r
library(Matrix)
library(dplyr)

# Use normalized expression matrix (log-normalized counts)
expr   <- GetAssayData(obj, layer = "data")  # genes x cells (use slot= for Seurat v4)
bc_vec <- obj$probe_barcode

# Filter to expressed genes (>= 5% of cells)
pct_expr  <- Matrix::rowMeans(expr > 0)
expr_filt <- expr[pct_expr >= 0.05, ]

groups <- sort(unique(bc_vec))
K  <- length(groups)
N  <- ncol(expr_filt)
df1 <- K - 1
df2 <- N - K

grand_mean <- Matrix::rowMeans(expr_filt)
SS_total   <- as.numeric(Matrix::rowSums(expr_filt^2) - N * grand_mean^2)
SS_between <- rep(0.0, nrow(expr_filt))

for (g in groups) {
  idx <- which(bc_vec == g)
  if (length(idx) < 2) next
  gm         <- as.numeric(Matrix::rowMeans(expr_filt[, idx, drop = FALSE]))
  SS_between <- SS_between + length(idx) * gm^2
}
SS_between <- SS_between - N * as.numeric(grand_mean)^2

eta_sq <- SS_between / SS_total
eta_sq[SS_total < 1e-10] <- 0
eta_sq[eta_sq < 0] <- 0
eta_sq[eta_sq > 1] <- 1

# F-statistic and p-value
SS_within <- pmax(SS_total - SS_between, 0)
F_stat    <- (SS_between / df1) / (SS_within / df2)
F_stat[SS_total < 1e-10] <- 0
p_val <- pf(F_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
p_adj <- p.adjust(p_val, method = "BH")

df_eta <- data.frame(
  gene          = rownames(expr_filt),
  eta_sq_bc     = eta_sq,
  F_stat        = F_stat,
  p_val         = p_val,
  p_adj         = p_adj,
  pct_expressed = as.numeric(pct_expr[pct_expr >= 0.05]),
  stringsAsFactors = FALSE
) %>% arrange(desc(eta_sq_bc)) %>% mutate(rank_bc = seq_len(n()))

cat(sprintf("Genes with eta^2 >= 1%%: %d\n", sum(df_eta$eta_sq_bc >= 0.01)))
cat(sprintf("Genes with eta^2 >= 5%%: %d\n", sum(df_eta$eta_sq_bc >= 0.05)))
print(head(df_eta[, c("gene", "eta_sq_bc", "rank_bc")], 10))
```

**Important caveat:** η² correlates with expression level due to statistical sensitivity — higher-expressed genes appear to have higher η² even if the true fold change is the same or smaller. The paper shows absolute log2FC *decreases* with expression (Supp. Fig. 2C), meaning the η²-expression correlation is a statistical artifact, not evidence of differential hybridization efficiency.

---

## FDP Analysis

Requires orthogonal barcoding (CRISPR guide assignments or hashtag oligos — see `experimental-design.md`).

```r
# --- Step 1: Ground truth true positives ---
# Compare biological conditions WITHIN a single probe barcode
bc_cells <- WhichCells(obj, expression = probe_barcode == "BC002")
obj_bc   <- subset(obj, cells = bc_cells)
Idents(obj_bc) <- obj_bc$guide_identity   # or hto_identity

de_within <- FindAllMarkers(
  obj_bc, only.pos = FALSE, min.pct = 0.1,
  logfc.threshold = 0, test.use = "wilcox"
)
true_positives <- unique(de_within$gene[
  !is.na(de_within$p_val_adj) &
  de_within$p_val_adj < 0.05 &
  abs(de_within$avg_log2FC) >= 0.3
])
cat(sprintf("True positive gene set: %d genes\n", length(true_positives)))

# --- Step 2: Confounded DE (cross-barcode, same biological condition) ---
guide_cells <- WhichCells(obj, expression = guide_identity == "guide_A")
obj_guide   <- subset(obj, cells = guide_cells)
Idents(obj_guide) <- obj_guide$probe_barcode

de_cross <- FindAllMarkers(
  obj_guide, only.pos = FALSE, min.pct = 0.1,
  logfc.threshold = 0, test.use = "wilcox"
)
sig_cross <- de_cross[!is.na(de_cross$p_val_adj) &
                      de_cross$p_val_adj < 0.05 &
                      abs(de_cross$avg_log2FC) >= 0.3, ]

# --- Step 3: FDP calculation ---
sig_cross$is_tp <- sig_cross$gene %in% true_positives
fdp_overall <- mean(!sig_cross$is_tp)
cat(sprintf("Overall FDP: %.1f%%\n", fdp_overall * 100))

# FDP at different log2FC thresholds (FDP sweep)
thresholds <- seq(0.1, 2.0, by = 0.1)
fdp_sweep <- sapply(thresholds, function(thresh) {
  sub_sig <- sig_cross[abs(sig_cross$avg_log2FC) >= thresh, ]
  if (nrow(sub_sig) == 0) return(NA)
  mean(!sub_sig$is_tp)
})
fdp_df <- data.frame(log2fc_thresh = thresholds, fdp = fdp_sweep)
```

---

## Visualization

### Volcano Plot

```r
library(ggplot2)
library(ggrepel)

FC_THRESH  <- 0.3
PVAL_THRESH <- 0.05

de_pair$sig <- with(de_pair,
  ifelse(p_val_adj < PVAL_THRESH & avg_log2FC >= FC_THRESH,  "up",
  ifelse(p_val_adj < PVAL_THRESH & avg_log2FC <= -FC_THRESH, "down", "ns")))

top_genes <- de_pair[de_pair$sig != "ns", ]
top_genes <- top_genes[order(abs(top_genes$avg_log2FC), decreasing = TRUE), ]
label_genes <- head(top_genes$gene, 15)

ggplot(de_pair, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.8, size = 1, shape = 16) +
  geom_vline(xintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed") +
  geom_hline(yintercept = -log10(PVAL_THRESH), linetype = "dashed") +
  scale_color_manual(values = c("up" = "tomato3", "down" = "steelblue3", "ns" = "grey")) +
  geom_text_repel(
    data = de_pair[de_pair$gene %in% label_genes, ],
    aes(label = gene), size = 3, color = "black", max.overlaps = 15
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Average log2FC", y = "-log10(p_val_adj)")
```

### Violin Plots of Top DE Genes

```r
# Top 4 genes by absolute log2FC
top4 <- head(de_pair[order(abs(de_pair$avg_log2FC), decreasing = TRUE), "gene"], 4)
VlnPlot(obj, features = top4, group.by = "probe_barcode", ncol = 2, pt.size = 0)
```

### Stacked Barplot of DE Gene Counts Per Barcode

```r
library(dplyr)
library(ggplot2)

deg_counts <- de_results[!is.na(de_results$p_val_adj), ] %>%
  mutate(direction = case_when(
    p_val_adj < 0.05 & avg_log2FC >= 0.3  ~ "Up",
    p_val_adj < 0.05 & avg_log2FC <= -0.3 ~ "Down",
    TRUE ~ "NS"
  )) %>%
  filter(direction != "NS") %>%
  count(cluster, direction)

ggplot(deg_counts, aes(x = cluster, y = n, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "tomato3", "Down" = "steelblue3")) +
  theme_classic() +
  labs(x = "Probe barcode", y = "Number of DE genes", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### FDP Sweep Plot

```r
ggplot(fdp_df[!is.na(fdp_df$fdp), ],
       aes(x = log2fc_thresh, y = fdp * 100)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50") +
  theme_classic() +
  labs(x = "|log2FC| threshold", y = "False Discovery Proportion (%)")
```

---

## Cross-Dataset Reproducibility Check

```r
# Load VCC reference
ref_de <- read.csv("supplementary_tables/Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv")
ref_bc <- ref_de[ref_de$cluster %in% c("BC002", "BC011"), ]
ref_avg <- aggregate(avg_log2FC ~ gene, data = ref_bc, FUN = mean)

# Merge with your pairwise DE
common <- merge(
  de_pair[, c("gene", "avg_log2FC")],
  ref_avg,
  by = "gene", suffixes = c("_yours", "_ref")
)

cor_result <- cor.test(common$avg_log2FC_yours, common$avg_log2FC_ref,
                       method = "spearman")
cat(sprintf("Spearman r = %.3f (p = %.2e)\n",
            cor_result$estimate, cor_result$p.value))
```

Paper benchmark: VCC vs PBMC BC002 vs BC011 Spearman ρ = **0.73**, p < 0.01 — showing the artifact is reproducible across completely independent datasets.

---

## Probe-Level h5 Analysis (Advanced — Python)

For investigating probe-level mechanisms using the CellRanger molecule_info.h5 output. Useful for understanding whether DE genes are driven by one probe vs all probes.

```python
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

CHUNK_SIZE = 50_000_000  # avoid memory overload

h5_path = "path/to/sample_molecule_info.h5"

with h5py.File(h5_path, "r") as f:
    # IMPORTANT: use probes/probe_id (not features/id)
    # probes/probe_id = full probe IDs like "ENSG00000000003|TSPAN6|8eab823"
    # features/id = Ensembl gene IDs (fewer rows, no probe-level info)
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

**Key finding from paper:** For the top DE gene SRP9, 2 probes exist — one is essentially silent (0.4% detection, 0.004 mean UMIs) and one is fully active (~22 UMIs/cell). The entire barcode DE effect on SRP9 is driven by the single active probe showing reduced expression in BC002. This probe-level heterogeneity points toward manufacturing batch effects in probe concentration rather than a systematic hybridization efficiency difference.

**Parse probe IDs** (format: `ENSG_ID-GENE_SYMBOL-PROBE_HASH`):
```python
# Extract gene symbol from probe ID
# Separator is "-" not "|"
import re
def parse_probe_id(probe_id):
    parts = probe_id.split("-")
    if len(parts) >= 2:
        return parts[1]  # gene symbol
    return probe_id

gene_symbols = [parse_probe_id(p) for p in probe_ids]
```
