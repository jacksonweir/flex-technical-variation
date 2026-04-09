---
name: flex-artifact-skill
description: Assesses probe barcode-associated technical variation in 10x Genomics Flex (Fixed RNA Profiling) scRNA-seq data. Use when user has Flex data with multiple probe barcodes and wants to check whether DE results are inflated by technical artifacts, or asks about probe barcode effects, sample barcoding variation, Flex v1 vs v2 technical differences, or spurious differential expression between barcodes. Based on Weir, Krebs, Chen 2026 (github.com/jacksonweir/flex-technical-variation).
metadata:
  author: Jackson Weir
  version: 1.0.0
---

# Probe Barcode Variation Assessment

## Overview

In 10x Genomics Flex v1, probe barcode identity drives substantial spurious differential expression (DE) between biologically identical cells — hundreds of genes, reproducible across lanes and datasets. This creates a high false discovery rate when probe barcodes are confounded with biological sample identity. Flex v2 significantly reduces this artifact, but the magnitude depends on whether probe hybridization was performed in bulk or independently per sample.

**Reference paper:** "Sample barcoding-associated technical variation in probe-based single-cell RNA sequencing" (Weir, Krebs, Chen 2026)
**Repository with analysis code and supplementary tables:** https://github.com/jacksonweir/flex-technical-variation

---

## Step 1 — Identify Assay Version and Experimental Design

**Flex v1**: Probe barcode identity is embedded in the probe hybridization chemistry. Barcodes appear as values like `BC001`–`BC016` in cell metadata. Technical variation between barcodes is inherent. **Any DE comparing conditions assigned to different barcodes will be confounded.**

**Flex v2**: Probe hybridization is barcode-free; barcoding happens in a separate oligo step. Barcodes follow a plate-well format (e.g. `A-A01`, `B-C04`) with the plate prefix varying depending on the plate used. The degree of technical variation depends on whether probe hybridization was bulk or per-sample.

Check the barcode naming in cell metadata to identify which assay was used. For detailed scenarios (Flex v1, Flex v2 bulk hyb, Flex v2 independent hyb) and how to identify which applies to your data, read:
`references/experimental-design.md`

---

## Step 2 — QC Metrics Per Probe Barcode

Before any DE analysis, verify QC metrics (cell count, median UMI, median genes detected, mitochondrial fraction) are comparable across barcodes. Systematically lower UMI counts or higher mitochondrial reads in specific barcodes are red flags. In the paper's VCC dataset, QC metrics per barcode were similar (Supp. Fig. 1), confirming the DE signal is not a QC artifact.

---

## Step 3 — Differential Expression Between Probe Barcodes

Run FindAllMarkers (one-vs-all, Wilcoxon, `min.pct = 0.1`, `logfc.threshold = 0`) with probe barcode as the identity. Count significant DE genes per barcode at `p_val_adj < 0.05` and `|log2FC| ≥ 0.3`.

Cross-reference your significant DE genes against the paper's supplementary tables to assess overlap with known probe barcode artifacts. Substantial overlap with the VCC (Table 1) or PBMC (Table 2) reference tables is strong evidence of probe barcode confounding.

**Recurrently probe barcode-affected genes** (high η², reproducible across independent datasets): **SRP9, SPCS1, ISCU, RRM2, SLC2A3, ITSN1, PHF6, ADCY2**

For code, read: `references/analysis-code.md` → Sections: FindAllMarkers, Pairwise DE

---

## Step 4 — Quantify Variance Explained (Eta-Squared)

Eta-squared (η²) measures the fraction of total expression variance attributable to probe barcode identity via one-way ANOVA: η² = SS_between / SS_total. Computed on NTC cells only (biologically equivalent cells, so any variance explained is purely technical). Run separately for probe barcode and lane to compare their magnitudes.

**Threshold from paper: η² ≥ 1% = probe barcode-affected gene**

Paper benchmarks (VCC NTC cells, all 3 lanes, all 16 probe barcodes):
- **519 genes** with η²_bc ≥ 1%
- **27 genes** with η²_lane ≥ 1%

For the full η² computation code (sparse matrix algebra, no densification required), read:
`references/analysis-code.md` → Section: Eta-Squared Computation

---

## Step 5 — Visualize the Effect

Key visualizations from the paper:
- **Volcano plots** (Fig. 1B, 1C): lane effect vs probe barcode effect side-by-side with shared axes
- **Violin plots** (Fig. 1D): per-cell expression of top η² genes across barcodes, split by lane — subtle, systematic shifts (no discrete on/off) are characteristic
- **Stacked barplot** (Fig. 1E): DE gene counts (up/down) per barcode from FindAllMarkers

If any of the recurrently probe barcode-affected genes (SRP9, SPCS1, ISCU, etc.) appear in your DE results, plotting their expression as violins across barcodes is a quick visual confirmation of the artifact.

For full plotting code, read: `references/analysis-code.md` → Section: Visualization

---

## Step 6 — FDP Analysis (Requires Orthogonal Barcoding)

If your cells have an independent identity label beyond probe barcode — such as **CRISPR guide assignments** or **hashtag oligonucleotides (HTOs)** from cell hashing — you can directly measure the false discovery proportion (FDP): the fraction of DE genes that are probe barcode artifacts.

The paper's approach (script `05_compute_fdp_sweep.R`): for each CRISPR target gene, run ground truth DE (NTC vs target cells, both within `TRUTH_BC`) and confounded DE (target cells in `CONFOUNDED_BC` vs NTC cells in `TRUTH_BC`). FDP = false positives / (false positives + true positives) at each |log2FC| threshold.

**Paper benchmarks (BC005 vs BC008, average-effect pair):**
- |log2FC| ≥ 0.3: median FDP = **67%**
- |log2FC| ≥ 1.0: median FDP = **34%**
- |log2FC| ≥ 1.7: FDP → **0%**

Worst-case pair (BC002 vs BC011): median FDP = **91%** at |log2FC| ≥ 0.3.

For code, read: `references/analysis-code.md` → Section: FDP Analysis

---

## Step 7 — Interpret and Recommend

For the full interpretation table (risk by experimental design), mitigation strategies, and Flex v1 vs Flex v2 comparison, read:
`references/interpretation.md`

**Quick risk assessment:**

| Situation | Risk |
|-----------|------|
| Flex v1, each biological condition on a separate barcode | HIGH |
| Flex v1, balanced across barcodes | MEDIUM |
| Flex v2, independent probe hyb per sample | REDUCED but non-zero |
| Flex v2, bulk probe hyb then split and barcoded | MINIMAL |

**Minimal mitigation if confounded:**
- Use `latent.vars = "probe_barcode"` with `test.use = "LR"` in FindMarkers
- Apply |log2FC| ≥ 1.7 threshold (FDP approaches 0% in paper benchmarks)
- Redesign: use Flex v2 with bulk probe hybridization, or distribute conditions across barcodes

---

## Reference Files

| File | When to read |
|------|-------------|
| `references/experimental-design.md` | To understand Flex v1 vs v2 chemistry, the three experimental scenarios, and how to identify which applies to your data |
| `references/analysis-code.md` | For complete R code: η² computation, FDP analysis, pairwise DE, all visualizations; Python code for probe-level h5 analysis |
| `references/interpretation.md` | For detailed risk table, paper benchmark numbers, problematic gene lists, mitigation strategies, cross-dataset reproducibility |

---

## Supplementary Tables

The paper's FindAllMarkers results for three datasets are in `supplementary_tables/`:

| File | Dataset | Cells used |
|------|---------|------------|
| `Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv` | VCC Flex v1, Lane 1 | All Lane 1 cells (all CRISPR perturbations, no guide filter); 16 barcodes |
| `Supplementary_Table_2_PBMC_Flexv1_FindAllMarkers.csv` | 10x PBMC 320k Flex v1 | All PBMC cells; 16 barcodes; FindAllMarkers one-vs-all |
| `Supplementary_Table_3_A375_Flexv2_FindAllMarkers.csv` | A375 Flex v2 | All cells across all conditions (8 ligands + PBS); 16 barcodes (2 per condition) |

Columns: `p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster (probe barcode), gene`

**Note:** Table 2 is the cleanest reference for identifying pure probe barcode artifacts — PBMC cells carry no CRISPR perturbations, so any DE between barcodes in this dataset is unambiguously technical. Table 1 includes all CRISPR guide populations, so a gene appearing in Table 1 could in principle be a genuine guide effect rather than a barcode artifact; cross-check against Table 2 or the problematic gene list to confirm.
