# Skill: Assess Probe Barcode-Associated Technical Variation in 10x Genomics Flex scRNA-seq

This skill guides you through assessing probe barcode-associated technical variation in 10x Genomics Flex (Fixed RNA Profiling) single-cell RNA-seq data, based on findings from:

**"Sample barcoding-associated technical variation in probe-based single-cell RNA sequencing"**
(Weir, Krebs, Chen — Broad Institute / Harvard)

Repository with analysis code and supplementary tables: `https://github.com/jacksonweir/flex-technical-variation`

---

## Overview

In 10x Genomics Flex v1, probe barcode identity drives substantial spurious differential expression (DE) between biologically identical cells — hundreds of genes, reproducible across lanes and datasets. This creates a high false discovery rate when probe barcodes are confounded with biological sample identity. Flex v2 significantly reduces this artifact, but the degree depends on whether probe hybridization is done in bulk (one tube) or independently per sample.

**When to use this skill:**
- You have 10x Genomics Flex (Fixed RNA Profiling) data
- Multiple probe barcodes are present (Flex v1 BC001–BC016, or Flex v2 sample barcodes)
- You plan to or have already performed DE analysis between conditions assigned to different probe barcodes
- You want to know whether your DE results may be inflated by technical artifacts

---

## Step 1 — Characterize Your Dataset

### 1a. Identify assay version

**Flex v1**: Up to 16 barcoded probe sets (BC001–BC016). The probe barcode is embedded directly in the probe hybridization chemistry — different barcodes use chemically distinct probe sets. Look for metadata columns named `probe_barcode` with values like `BC002`, `BC011`, etc.

**Flex v2**: A barcode-free probe set hybridizes first; sample barcoding is done in a separate oligo hybridization step afterward. Look for metadata with values like `A-A01` through `A-H12` (96-plex format).

```r
# Check what probe barcode values look like
table(obj$probe_barcode)
# or check metadata column names
colnames(obj@meta.data)
```

### 1b. Understand your experimental design — CRITICAL

The magnitude of the technical artifact depends entirely on **when probe hybridization was done relative to sample barcoding**:

**Scenario A — Flex v1 (any design)**
Probe barcode is inseparable from probe hybridization. Technical variation between barcodes is inherent to the chemistry. Expect substantial spurious DE (see Step 3 for quantification). Any DE analysis comparing conditions assigned to different barcodes will be confounded.

**Scenario B — Flex v2, bulk probe hybridization (one tube), then split and barcoded**
Cells from all conditions are pooled and probe-hybridized together, then split and assigned to different barcodes. This decouples probe hybridization from barcoding and approaches zero technical variation between barcodes (see Supp. Fig. 5 in the paper). DE analysis is reliable.

**Scenario C — Flex v2, independent probe hybridization per sample**
Each sample/condition is probe-hybridized separately, then barcoded. Reduced variation compared to Flex v1, but non-zero. The A375 ligand dictionary dataset (Fig. 2B in the paper) is this scenario — still shows some spurious DE, especially for mitochondrial genes and other technically sensitive genes.

**Ask yourself or check your protocol:**
- Were all samples pooled before probe hybridization? → Scenario B (low risk)
- Was each sample hybridized separately? → Scenario A or C (proceed with caution)

### 1c. Check for orthogonal barcoding

Do your cells have an additional, independent identity label beyond the probe barcode? For example:
- **CRISPR guide assignments** (guide RNA identity)
- **Hashtag oligonucleotides (HTOs)** from CITE-seq/cell hashing
- **Known cell line identity** (mixed species or genomic variants)

If yes, you can perform a **false discovery proportion (FDP) analysis** to directly quantify what fraction of your DE results are artifacts (see Step 6).

---

## Step 2 — QC Metrics Per Probe Barcode

Before any DE analysis, verify that QC metrics are comparable across probe barcodes. Systematic QC differences could either mask or inflate the technical variation signal.

```r
library(Seurat)
library(ggplot2)
library(dplyr)

# Compute per-barcode QC summaries
qc_summary <- obj@meta.data %>%
  group_by(probe_barcode) %>%
  summarise(
    n_cells       = n(),
    median_nCount = median(nCount_RNA),
    median_nFeat  = median(nFeature_RNA),
    median_mito   = median(percent.mito),
    .groups = "drop"
  )
print(qc_summary)

# Violin plots of nCount, nFeature, percent.mito per barcode
VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"),
        group.by = "probe_barcode", ncol = 3, pt.size = 0)
```

**What to look for:**
- Cell counts should be reasonably balanced across barcodes (within ~2–3×)
- Median UMI counts and gene counts should be similar across barcodes
- Percent mitochondrial reads should be similar; systematic differences in mito percentage are a red flag (seen in Supp. Fig. 4 for Flex v2 A375 dataset)

In the paper's VCC dataset, QC metrics per barcode were similar (Supp. Fig. 1), confirming the DE signal is not a QC artifact.

---

## Step 3 — Differential Expression Between Probe Barcodes

### 3a. Run FindAllMarkers across all probe barcodes

This tests each probe barcode against all others pooled (one-vs-all). Parameters match those used in the paper:

```r
Idents(obj) <- obj$probe_barcode

de_results <- FindAllMarkers(
  obj,
  only.pos        = FALSE,
  min.pct         = 0.1,
  logfc.threshold = 0,
  test.use        = "wilcox",
  verbose         = TRUE
)

# Count significant DE genes per barcode
sig <- de_results[!is.na(de_results$p_val_adj) &
                  de_results$p_val_adj < 0.05 &
                  abs(de_results$avg_log2FC) >= 0.3, ]
print(sort(table(sig$cluster)))
```

### 3b. Run pairwise DE for a specific barcode pair (optional)

For a focused comparison, especially if you suspect a particular barcode pair is problematic:

```r
# Example: compare BC002 vs BC011 (worst-case pair in the paper)
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

### 3c. Cross-reference with the paper's supplementary tables

The paper provides FindAllMarkers results for three datasets in `supplementary_tables/`:

| File | Dataset | Assay | Notes |
|------|---------|-------|-------|
| `Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv` | VCC (Arc Institute) | Flex v1 | NTC cells only, Lane 1; 16 barcodes |
| `Supplementary_Table_2_PBMC_Flexv1_FindAllMarkers.csv` | 10x PBMC 320k | Flex v1 | 16 barcodes |
| `Supplementary_Table_3_A375_Flexv2_FindAllMarkers.csv` | A375 ligand dictionary | Flex v2 | 8 barcodes, biologically distinct conditions |

```r
# Load the Flex v1 reference (VCC dataset)
vcc_de <- read.csv("supplementary_tables/Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv")
# Columns: p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster (probe barcode), gene

# Significant DE genes per barcode in the reference dataset
vcc_sig <- vcc_de[!is.na(vcc_de$p_val_adj) &
                  vcc_de$p_val_adj < 0.05 &
                  abs(vcc_de$avg_log2FC) >= 0.3, ]
sort(table(vcc_sig$cluster))

# Find genes that appear in BOTH your dataset and the reference
your_sig_genes <- unique(sig$gene)  # from your FindAllMarkers
ref_sig_genes  <- unique(vcc_sig$gene)
overlap <- intersect(your_sig_genes, ref_sig_genes)
cat(sprintf("Genes in your DE results that are also in the VCC reference: %d\n",
            length(overlap)))
```

**Interpretation:**
- If your DE gene list overlaps substantially with the paper's VCC or PBMC supplementary tables, those genes are likely probe barcode artifacts — they show up in biologically unrelated datasets.
- The paper found 217 significant DE genes between BC002 and BC011 in Lane 1 of the VCC dataset using all cells (|log2FC| > 0.3, p_adj < 0.05).
- Key recurrently affected genes: **SRP9, SPCS1, ISCU, RRM2, SLC2A3** — if these appear in your DE results, they are almost certainly technical.

---

## Step 4 — Visualize the Effect

### 4a. Volcano plots per barcode pair

```r
library(ggplot2)
library(ggrepel)

FC_THRESH <- 0.3
PVAL_THRESH <- 0.05

# Use de_pair from Step 3b
de_pair$sig <- with(de_pair,
  ifelse(p_val_adj < PVAL_THRESH & avg_log2FC >= FC_THRESH, "up",
  ifelse(p_val_adj < PVAL_THRESH & avg_log2FC <= -FC_THRESH, "down", "ns")))

# Genes to label
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

### 4b. Violin plots of top DE genes

```r
# Top 4 genes by absolute log2FC
top4 <- head(de_pair[order(abs(de_pair$avg_log2FC), decreasing = TRUE), "gene"], 4)

VlnPlot(obj, features = top4, group.by = "probe_barcode",
        ncol = 2, pt.size = 0)
```

**What to look for in violin plots:**
- Do the expression differences look biologically meaningful (discrete on/off, large fold changes) or are they subtle shifts in distribution across all barcodes?
- Subtle, systematic shifts across many barcodes with no clear on/off pattern are characteristic of probe barcode technical variation (see Fig. 1D in the paper for SRP9, SPCS1, ISCU).

### 4c. Stacked barplot of DE gene counts

```r
# Counts of up/down DEGs per barcode
deg_counts <- de_results[!is.na(de_results$p_val_adj), ] %>%
  mutate(
    direction = case_when(
      p_val_adj < 0.05 & avg_log2FC >= 0.3 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC <= -0.3 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(direction != "NS") %>%
  count(cluster, direction)

ggplot(deg_counts, aes(x = cluster, y = n, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "tomato3", "Down" = "steelblue3")) +
  theme_classic() +
  labs(x = "Probe barcode", y = "Number of DE genes", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

---

## Step 5 — Quantify Variance Explained (Eta-Squared)

Eta-squared (η²) measures the fraction of total expression variance attributable to probe barcode identity. Computed via one-way ANOVA: η² = SS_between / SS_total.

**Interpretation threshold from the paper: η² ≥ 1% = probe barcode-affected gene.**

In the VCC Flex v1 dataset: 519 genes reached η² ≥ 1% for probe barcode vs only 27 genes for lane identity — showing the effect is barcode-specific, not a general batch effect.

```r
library(Matrix)

# Use the normalized expression matrix (log-normalized counts)
expr    <- GetAssayData(obj, layer = "data")  # genes x cells; use slot= for Seurat v4
bc_vec  <- obj$probe_barcode

# Filter to expressed genes (>= 5% of cells)
pct_expr   <- Matrix::rowMeans(expr > 0)
expr_filt  <- expr[pct_expr >= 0.05, ]

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
)
df_eta <- df_eta[order(df_eta$eta_sq_bc, decreasing = TRUE), ]

cat(sprintf("Genes with eta^2 >= 1%%: %d\n", sum(df_eta$eta_sq_bc >= 0.01)))
cat(sprintf("Genes with eta^2 >= 5%%: %d\n", sum(df_eta$eta_sq_bc >= 0.05)))
print(head(df_eta[, c("gene", "eta_sq_bc")], 10))
```

**Important caveat:** η² is correlated with expression level — higher-expressed genes are more statistically sensitive and will appear to have higher η² even if the true effect size (fold change) is the same. The paper shows that absolute log2FC *decreases* with expression (Supp. Fig. 2C), meaning the η²-expression correlation is a statistical artifact, not evidence that probe hybridization efficiency drives the effect.

---

## Step 6 — FDP Analysis (Requires Orthogonal Barcoding)

If your cells have an independent identity label beyond the probe barcode (CRISPR guide, hashtag/HTO, known cell line), you can directly measure the **false discovery proportion (FDP)** — the fraction of your DE genes that are probe barcode artifacts rather than true biology.

**Concept:**
- **Within-barcode DE** (same probe barcode, different biological conditions) = ground truth true positives
- **Cross-barcode DE** (same biological condition, different probe barcodes) = ground truth false positives
- For each DE gene in your actual confounded analysis, label it TP if it appears in the within-barcode analysis, FP if it only appears in the cross-barcode analysis

```r
# Assumes: obj has both 'probe_barcode' (e.g. BC002) and 'guide_identity' (e.g. guide_A) metadata

# --- Within-barcode DE (ground truth) ---
# Pick one barcode, compare biological conditions within it
bc_cells <- WhichCells(obj, expression = probe_barcode == "BC002")
obj_bc   <- subset(obj, cells = bc_cells)
Idents(obj_bc) <- obj_bc$guide_identity

de_within <- FindAllMarkers(
  obj_bc, only.pos = FALSE, min.pct = 0.1,
  logfc.threshold = 0, test.use = "wilcox"
)
true_positives <- unique(de_within$gene[
  !is.na(de_within$p_val_adj) &
  de_within$p_val_adj < 0.05 &
  abs(de_within$avg_log2FC) >= 0.3
])

# --- Cross-barcode DE (confounded analysis) ---
# Compare same biological condition across different probe barcodes
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

# --- FDP calculation ---
sig_cross$is_tp <- sig_cross$gene %in% true_positives
fdp_overall <- mean(!sig_cross$is_tp)
cat(sprintf("FDP: %.1f%% of DE genes are probe barcode artifacts\n", fdp_overall * 100))
```

**Paper benchmarks** (BC005 vs BC008, average-effect pair; NTC cells, VCC dataset):
- At |log2FC| ≥ 0.3: median FDP = **67%**
- At |log2FC| ≥ 1.0: median FDP = **34%**
- At |log2FC| ≥ 1.7: FDP drops to **0%**

For the worst-case pair (BC002 vs BC011): median FDP = **91%** at |log2FC| ≥ 0.3.

---

## Step 7 — Cross-Dataset Reproducibility

A key feature of the probe barcode effect is its **reproducibility across independent datasets**. Genes affected in one dataset tend to be affected in another using the same probe barcodes.

```r
# Load reference DE results from this repository
ref_de <- read.csv("supplementary_tables/Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv")
ref_de_bc002_bc011 <- ref_de[ref_de$cluster %in% c("BC002", "BC011"), ]

# Your pairwise DE (from Step 3b, same barcode pair)
# de_pair from FindMarkers, gene column added

# Merge and compare log2FC
common <- merge(
  de_pair[, c("gene", "avg_log2FC", "p_val_adj")],
  # Approximate: take average log2FC for BC002 vs BC011 from ref
  aggregate(avg_log2FC ~ gene, data = ref_de_bc002_bc011, FUN = mean),
  by = "gene", suffixes = c("_yours", "_ref")
)

cor_result <- cor.test(common$avg_log2FC_yours, common$avg_log2FC_ref,
                       method = "spearman")
cat(sprintf("Spearman r between your log2FC and reference: %.3f (p=%.2e)\n",
            cor_result$estimate, cor_result$p.value))
```

**What to look for:** A significant positive correlation means the DE pattern in your data mirrors the probe barcode technical artifact seen in the reference dataset — strong evidence that your DE results are technically confounded.

---

## Step 8 — Interpretation and Recommendations

### If you have Flex v1 data:

**The probe barcode effect is present in your data.** The question is whether it confounds your biological comparisons.

| Your experimental design | Risk level | Recommendation |
|--------------------------|------------|----------------|
| Biological conditions each assigned to separate probe barcodes | **HIGH** | Do not trust DE between conditions without mitigation |
| Biological conditions distributed across multiple barcodes (balanced) | **MEDIUM** | Include probe barcode as a covariate in regression-based DE |
| Same biological condition across all barcodes (NTC/superloading) | **Low (for that comparison)** | Probe barcode DE is pure artifact; use as negative control |

**Mitigation strategies (from the paper):**
1. Use `FindMarkers` with `latent.vars = "probe_barcode"` to regress out the barcode effect (requires LR or mixed model test)
2. Use pseudobulk methods (DESeq2, edgeR) with probe barcode as a blocking factor
3. Avoid |log2FC| < 1.5–2.0 as threshold when barcodes are confounded — the paper shows FDP approaches 0% at |log2FC| ≥ 1.7

### If you have Flex v2 data with bulk probe hybridization:

Minimal risk. The technical variation is near zero when probe hybridization is done in bulk before splitting. Proceed normally.

### If you have Flex v2 data with independent probe hybridization:

Reduced risk compared to Flex v1, but non-zero. Check QC metrics (Step 2) carefully — mitochondrial gene DE between conditions is a red flag. The A375 dataset in this paper (Supp. Fig. 4) shows that even with independent probe hyb, the effect is substantially reduced vs Flex v1. Use the A375 supplementary table (`Supplementary_Table_3`) as a reference for which genes are artifactually DE in Flex v2.

### Key genes to flag as likely technical artifacts:

These genes appear in the paper's VCC and/or PBMC Flex v1 reference datasets and recurrently show probe barcode-associated DE:

**Flex v1 high-risk genes** (η² ≥ 5% in VCC dataset): SRP9, SPCS1, ISCU, RRM2, SLC2A3, ITSN1, PHF6, ADCY2

**Flex v2 flag genes** (top DE in A375 independent hyb, not biologically driven): Check Table 3 for genes appearing in multiple conditions with no clear biological rationale — mitochondrial genes (MT-ND1, MT-CO1, etc.) are particularly prominent in Flex v2 Supp. Fig. 4.

---

## Step 9 — Probe-Level Analysis (Advanced)

The probe barcode effect operates primarily at the probe level. For multi-probe genes, the DE signal is often driven by a single probe showing dropout or reduced efficiency in specific barcodes — while sibling probes of the same gene may be unaffected (Fig. 2E–F in the paper).

If you have access to the probe-level molecule_info.h5 output from CellRanger:

```python
# Python — requires h5py, numpy, pandas, scipy
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

CHUNK_SIZE = 50_000_000  # read in chunks to avoid memory overload

h5_path = "path/to/sample_molecule_info.h5"

with h5py.File(h5_path, "r") as f:
    # Probe IDs are in probes/probe_id (not features/id)
    probe_ids   = f["probes/probe_id"][:].astype(str)
    barcodes_h5 = f["barcodes"][:].astype(str)

    n_total = f["count"].shape[0]
    probe_idx_all   = []
    barcode_idx_all = []
    count_all       = []

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

**What to look for:**
- For a DE gene (e.g., SRP9), extract per-probe, per-barcode expression and check whether all probes show the same pattern or only one probe is responsible.
- In the paper, SRP9 has 2 probes — one is essentially silent (0.4% detection), the other is fully active (~22 UMIs/cell). The entire barcode DE effect on SRP9 is driven by the single active probe showing reduced expression in BC002.
- This probe-level heterogeneity suggests the mechanism involves probe-level concentration or synthesis differences during manufacturing, not a systematic hybridization efficiency difference.

---

## Quick-Start Checklist

```
[ ] 1. Identify assay version (Flex v1 / v2) and probe barcode metadata column
[ ] 2. Understand experimental design: was probe hyb done in bulk or per-sample?
[ ] 3. Check QC metrics per barcode — any systematic differences?
[ ] 4. Run FindAllMarkers across probe barcodes (Step 3a)
[ ] 5. Count significant DE genes per barcode; flag barcodes with high counts
[ ] 6. Cross-reference your top DE genes with supplementary_tables/ in this repo
[ ] 7. Check if key sentinel genes (SRP9, SPCS1, ISCU) appear in your results
[ ] 8. If orthogonal barcoding available: run FDP analysis (Step 6)
[ ] 9. Compute eta-squared to quantify which genes are most affected (Step 5)
[ ] 10. Apply appropriate mitigation based on your experimental design (Step 8)
```

---

## Reference

For full methods, statistics, and figures, see the paper and the analysis code in this repository. The R scripts in `scripts/` are numbered and annotated; follow the ordering in the README to reproduce all figures.

Key scripts for the analyses described above:
- `scripts/figure1/compute/02_compute_pairwise_de.R` — pairwise DE (BC002 vs BC011)
- `scripts/figure1/compute/03_compute_findallmarkers_lane1.R` — FindAllMarkers (Flex v1)
- `scripts/figure1/compute/04_compute_variance_explained.R` — η² computation
- `scripts/figure1/compute/05_compute_fdp_sweep.R` — FDP sweep analysis
- `scripts/figure2/compute/09_compute_a375_de.R` — FindAllMarkers (Flex v2)
- `scripts/figure2/compute/10_compute_variance_explained_a375.R` — η² for Flex v2
