# Interpretation and Recommendations Reference

## Paper Benchmark Numbers

### Flex v1 — VCC Dataset (NTC cells, all 16 barcodes, 3 lanes)

| Metric | Value |
|--------|-------|
| DE genes BC002 vs BC011 (worst-case pair, all cells) | **217** at \|log2FC\| > 0.3 |
| Genes with η² ≥ 1% (probe barcode) | **519** |
| Genes with η² ≥ 1% (lane identity) | **27** |
| FDP at \|log2FC\| ≥ 0.3 (BC002 vs BC011) | **~91%** |
| FDP at \|log2FC\| ≥ 0.3 (BC005 vs BC008, average pair) | **67%** |
| FDP at \|log2FC\| ≥ 1.0 (BC005 vs BC008) | **34%** |
| FDP → 0% at | \|log2FC\| ≥ **1.7** |
| Cross-dataset reproducibility (VCC vs PBMC, Spearman ρ) | **0.73** (p < 0.01) |

### Flex v2 — A375 Dataset (8 ligand conditions, independent probe hyb)

- Substantially fewer DE genes per barcode than Flex v1 (see Fig. 2B vs Fig. 1E)
- Some residual DE, especially for mitochondrial genes and housekeeping genes
- Sentinel genes most affected in Flex v2: check `supplementary_tables/Supplementary_Table_3_A375_Flexv2_FindAllMarkers.csv`

### Flex v2 — Bulk Probe Hybridization (Supp. Fig. 5)

- DE between barcodes approaches zero when cells are probe-hybridized in a single pooled tube
- This is the minimal-artifact scenario

---

## Risk Table by Experimental Design

| Experimental design | Assay | Risk | Action |
|--------------------|-------|------|--------|
| Each biological condition on a separate probe barcode | Flex v1 | **HIGH** | Do not trust DE without mitigation; median FDP 67–91% at \|log2FC\| 0.3 |
| Biological conditions split across multiple barcodes (balanced) | Flex v1 | **MEDIUM** | Include probe barcode as covariate; use pseudobulk methods |
| Same biological condition in all barcodes (NTC/superloading) | Flex v1 | N/A — pure artifact | Use as negative control to characterize the artifact |
| Each condition on separate barcode, independent probe hyb | Flex v2 | **REDUCED** | Check empirically; mitochondrial and housekeeping genes most at risk |
| Cells pooled before probe hyb, then split and barcoded | Flex v2 | **MINIMAL** | Proceed normally; run a quick DE check to confirm |

---

## Sentinel Artifact Genes

These genes appear recurrently in probe barcode DE results across independent datasets. If any of these appear in your DE gene list, they are almost certainly technical artifacts:

**Flex v1 high-η² genes** (η² ≥ 5% in VCC dataset):
`SRP9, SPCS1, ISCU, RRM2, SLC2A3, ITSN1, PHF6, ADCY2`

**Additional Flex v1 affected genes** (η² ≥ 1%): 519 total — check `supplementary_tables/Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv` for the full list.

**Functional categories enriched among affected genes:**
- RNA processing and splicing (snRNP complexes)
- Mitochondrial respiration and inner membrane
- Translation regulation
- (Based on GO ORA of 519 η²-affected genes)

---

## Mitigation Strategies

### 1. Regress out probe barcode in DE (within Seurat)

```r
# Requires LR test (or mixed model); not available with wilcox
de_corrected <- FindMarkers(
  obj,
  ident.1        = "condition_A",
  ident.2        = "condition_B",
  test.use       = "LR",
  latent.vars    = "probe_barcode",
  min.pct        = 0.1,
  logfc.threshold = 0
)
```

**Note:** The paper found that UMI correction via LR test does not fully eliminate the artifact (Spearman ρ = 0.82 between corrected and uncorrected log2FC for affected genes). Use with caution.

### 2. Pseudobulk methods with blocking factor

```r
library(DESeq2)

# Aggregate counts per (sample, barcode) combination
# then use barcode as blocking factor in the model

# Create pseudobulk counts
pseudo_counts <- AggregateExpression(obj,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = c("sample_id", "probe_barcode")
)$RNA

# DESeq2 with probe_barcode as blocking factor
dds <- DESeqDataSetFromMatrix(
  countData = pseudo_counts,
  colData   = sample_metadata,
  design    = ~ probe_barcode + condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "A", "B"))
```

### 3. Raise the log2FC threshold

Based on FDP sweep analysis: the artifact contributes minimally above |log2FC| = 1.7 in Flex v1.

- |log2FC| ≥ 0.3: median FDP ~67% (average barcode pair)
- |log2FC| ≥ 1.0: median FDP ~34%
- |log2FC| ≥ 1.7: FDP → 0%

This is a conservative but effective filter for exploratory analyses.

### 4. Experimental design correction (preferred)

If possible, redesign the experiment:
- Use Flex v2 with bulk probe hybridization (Scenario B) — minimal artifact
- Distribute biological conditions across multiple barcodes (not one condition per barcode)
- Include NTC cells in each barcode as an internal negative control

---

## Interpreting Cross-Dataset Overlap

When your DE genes substantially overlap with the paper's reference datasets (VCC Flex v1 or PBMC Flex v1), this is strong evidence of technical confounding:

```
Your DE genes ∩ VCC reference → likely probe barcode artifact
Your DE genes ∩ A375 Flex v2 reference → may persist in Flex v2
Your DE genes ∩ neither reference → could be biology (still check manually)
```

The cross-dataset reproducibility (Spearman ρ = 0.73 for VCC vs PBMC BC002/BC011 fold changes) shows the artifact is a consistent feature of specific probe barcode pairs, not random noise.

---

## Mechanistic Notes (What We Know and Don't Know)

**What drives the effect (current understanding):**
- Probe barcode technical variation is not explained by probe sequence features (GC content does not predict η² after controlling for expression level)
- Synthesis sequence errors in probe barcodes are NOT the mechanism (BAM mismatch rates are uncorrelated with expression ranks across barcodes)
- The most likely mechanism is batch-level concentration or quality differences in probe manufacturing — some barcodes have slightly lower probe concentrations or higher truncation rates, leading to globally reduced capture efficiency
- The effect is global (affects all genes proportionally to their expression level), not gene-specific

**Probe-level heterogeneity:**
- For multi-probe genes, DE is often driven by a single probe rather than all probes of that gene
- This is consistent with probe-level manufacturing batch effects rather than a gene-level mechanism
- For Flex v1 VCC: 561/1,366 target probes (41%) have η² ≥ 1%; 469/507 genes (93%) have at least 1 affected probe

**Flex v2 improvement:**
- Flex v2 decouples probe hybridization from barcoding, eliminating the probe barcode identity as a source of probe set variation
- The residual effect in Flex v2 with independent probe hyb is likely due to tube-to-tube variability in the separate barcoding step, which is smaller in magnitude than Flex v1's intrinsic probe-barcode coupling
