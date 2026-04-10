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

> **Not all steps are necessary for every dataset.** Start with the highest-yield steps:
> identify the experimental design (Step 1), check QC metrics (Step 2), and run the VCC
> reference comparison (Steps 3–3b). These three steps together will usually give a clear
> answer about whether the artifact is present and how severe it is. Ask the user whether
> to proceed to eta-squared computation (Step 4), visualization (Step 5), or FDP analysis
> (Step 6) based on what the initial steps reveal and what level of characterization is needed.

---

## Step 0 — Assign Probe Barcode Identity to Cells

Before any analysis, confirm that each cell in your Seurat object has a `probe_barcode`
metadata column (values like `BC001`–`BC016` for Flex v1). This column is required for
all downstream steps.

**If the column is already present** (e.g., loaded from per-sample cellranger multi outputs
and labeled during object construction), proceed to Step 1.

**If the column is absent** (e.g., you loaded a combined or aggregated matrix without
per-sample metadata), extract probe barcode identity directly from the cell barcode sequences.
In Flex v1, the last 8 bp of the 16 bp cell barcode encode the probe barcode identity. Match
these 8 bp sequences against the BC001–BC016 whitelist to annotate each cell.

Full code and the sequence whitelist are in: `references/extract-probe-barcodes.md`

This approach is implemented in `scripts/figure1/01_setup_vcc_seurat.R` in this repository.

---

## Step 1 — Identify Assay Version and Experimental Design

**Flex v1**: Probe barcode identity is embedded in the probe hybridization chemistry. Barcodes appear as values like `BC001`–`BC016` in cell metadata. Technical variation between barcodes is inherent. **Any DE comparing conditions assigned to different barcodes will be confounded.**

**Flex v2**: Probe hybridization is barcode-free; barcoding happens in a separate oligo step. Barcodes follow a plate-well format (e.g. `A-A01`, `B-C04`) with the plate prefix varying depending on the plate used. The degree of technical variation depends on whether probe hybridization was bulk or per-sample.

### Confounded Designs — Very Common and High Risk

**A critical and common situation**: in many Flex v1 experiments, each biological
condition is placed on a separate probe barcode. For example, treated cells on BC002 and
control cells on BC011, or patient A on BC001 and patient B on BC002. In this case,
**probe barcode identity is perfectly confounded with sample identity**. Any DE analysis
between your biological conditions will produce a mix of two signals:

1. Genuine biological differences between the conditions
2. Technical artifact genes introduced by the probe barcode itself

These two sources of DE are completely inseparable without an orthogonal label (e.g.,
CRISPR guide assignments or HTOs — see Step 6). The paper benchmarks this directly: at
|log2FC| ≥ 0.3, the median false discovery proportion is **67–91%** depending on which
barcode pair is used. This means the majority of "significant" DE genes in a confounded
comparison may not reflect biology at all.

**For confounded designs, the immediate priority is the VCC reference comparison (Step 3b):**
run FindAllMarkers on your data with probe barcode as the identity, then compare each
barcode's one-vs-all DE to the VCC reference. This directly tests whether the same artifact
genes are inflating your results.

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

## Step 3b — VCC Reference Comparison and Confounded Design Check

**This step is especially important when probe barcode is perfectly confounded with sample
identity** (e.g., all treated cells on BC002, all control cells on BC011). In that design,
probe barcode technical artifacts are indistinguishable from biological DE.

**Primary approach — one-vs-all per barcode:** Run FindAllMarkers on your data with probe
barcode as the identity. For each barcode, scatter-plot your avg_log2FC (barcode vs all
others) against the VCC reference avg_log2FC for the same barcode from Supplementary Table
1. Since both use the same one-vs-all design, the fold-changes are directly comparable. A
positive Spearman ρ for any barcode means that barcode's artifact profile is present in
your data. This is the most direct and interpretable diagnostic for confounded designs.

**Advanced — pairwise comparison:** For a specific barcode pair of interest (e.g., the two
barcodes that correspond to your treatment vs control), run pairwise DE and compare to an
approximated VCC pairwise log2FC derived from Table 1.

**Supplementary Table 1 (VCC) is the primary reference** because VCC NTC cells are
biologically identical across all 16 barcodes — any DE between barcodes is unambiguously
technical. Table 2 (PBMC) provides a second independent confirmation.

Both tables are in `references/` in this skill directory for direct access.

For full code (one-vs-all scatter plots, pairwise comparison, overlap statistics), read:
`references/vcc-comparison.md`

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
| `references/extract-probe-barcodes.md` | If your Seurat object lacks a probe barcode metadata column — extract BC001–BC016 identity from cell barcode sequences |
| `references/experimental-design.md` | To understand Flex v1 vs v2 chemistry, the three experimental scenarios, and how to identify which applies to your data |
| `references/analysis-code.md` | For complete R code: η² computation, FDP analysis, pairwise DE, all visualizations; Python code for probe-level h5 analysis |
| `references/interpretation.md` | For detailed risk table, paper benchmark numbers, problematic gene lists, mitigation strategies, cross-dataset reproducibility |
| `references/vcc-comparison.md` | For scatter plot code comparing your DE results to the VCC reference, and the confounded design check |

---

## Supplementary Tables

The paper's FindAllMarkers results for three datasets are in `supplementary_tables/`:

| File | Dataset | Cells used |
|------|---------|------------|
| `Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv` | VCC Flex v1, Lane 1 | All Lane 1 cells (all CRISPR perturbations, no guide filter); 16 barcodes |
| `Supplementary_Table_2_PBMC_Flexv1_FindAllMarkers.csv` | 10x PBMC 320k Flex v1 | All PBMC cells; 16 barcodes; FindAllMarkers one-vs-all |
| `Supplementary_Table_3_A375_Flexv2_FindAllMarkers.csv` | A375 Flex v2 | All cells across all conditions (8 ligands + PBS); 16 barcodes (2 per condition) |

Columns: `p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster (probe barcode), gene`

**Table 1 (VCC) is the primary reference for the confounded design check.** VCC NTC cells
are biologically identical across all 16 barcodes, making this the gold-standard artifact
gene list. Tables 1 and 2 are also available in `references/` within this skill directory
for direct use without navigating to `supplementary_tables/`.

**Table 2 (PBMC)** carries no CRISPR perturbations, so any DE between barcodes is
unambiguously technical. Use it as a second independent confirmation when a gene appears
in Table 1 — genes significant in both VCC and PBMC are highly reliable artifact indicators.

**Note on Table 1 and CRISPR guide effects:** Table 1 includes all CRISPR guide populations.
A gene appearing in Table 1 could in principle reflect a genuine CRISPR perturbation rather
than a barcode artifact. Cross-check against Table 2 or the recurrently affected gene list
(SRP9, SPCS1, ISCU, RRM2, SLC2A3, ITSN1, PHF6, ADCY2) to confirm technical origin.
