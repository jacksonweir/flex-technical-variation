# Experimental Design Reference

## Flex v1 vs Flex v2 Chemistry

### Flex v1
- Up to 16 barcoded probe sets (BC001–BC016)
- The probe barcode sequence is embedded **within the probe itself** — different barcodes use chemically distinct probe sets
- Probe hybridization and sample barcoding happen simultaneously in the same step
- Technical variation between barcodes is inherent to the chemistry and cannot be avoided
- Metadata values: `BC001`, `BC002`, ..., `BC016`

### Flex v2
- A **single barcode-free probe set** hybridizes to all samples in one step
- Sample barcoding is done in a **separate oligo hybridization step** afterward
- This decoupling is why Flex v2 shows substantially reduced barcode-associated technical variation
- Metadata values typically: `A-A01` through `A-H12` (96-plex format)

---

## Three Experimental Scenarios

### Scenario A — Flex v1 (any experimental design)

Probe barcode is inseparable from probe hybridization. Technical variation between barcodes is inherent.

**Risk level: HIGH for any biological comparison between barcodes.**

Even if the experiment is well-designed biologically, probe barcode technical effects will confound DE between conditions assigned to different barcodes.

**How to identify:** Metadata column contains `BC001`–`BC016`-style values.

---

### Scenario B — Flex v2, bulk probe hybridization

All cells from all conditions are **pooled into one tube**, probe-hybridized together, then split and assigned to different barcodes afterward.

- Since all cells see the same probe hybridization environment, probe-level technical variation is eliminated
- DE between barcodes should be near-zero technical artifact (see Supp. Fig. 5 in paper)
- This is the safest Flex v2 design for comparative experiments

**Risk level: MINIMAL.**

**How to identify:** Ask the experimentalist. Protocol notes will say cells were pooled before probe hybridization. Often used for superloading experiments where one cell type is split across barcodes.

---

### Scenario C — Flex v2, independent probe hybridization per sample

Each biological sample is probe-hybridized **separately**, then barcoded. This mirrors how Flex v1 is typically used.

- Reduced variation compared to Flex v1 (because probe sets are identical across samples)
- But non-zero: independent hyb introduces tube-to-tube variability in probe concentration, hybridization efficiency, etc.
- The A375 ligand dictionary dataset (Fig. 2B in paper) is this scenario — still shows some spurious DE, especially for mitochondrial and housekeeping genes
- Check `supplementary_tables/Supplementary_Table_3_A375_Flexv2_FindAllMarkers.csv` for which genes are most affected in Flex v2

**Risk level: REDUCED but non-zero. Assess empirically.**

**How to identify:** Ask the experimentalist. Separate probe hyb per condition, then pooled for barcoding.

---

## How to Identify Your Scenario

1. Check metadata column names and values:
   ```r
   colnames(obj@meta.data)
   table(obj$probe_barcode)
   ```

2. If values are `BC001`–`BC016`: you have **Flex v1** (Scenario A).

3. If values look like `A-A01`, `A-B01`, etc.: you have **Flex v2**.
   - Then ask: were cells pooled before probe hybridization (Scenario B) or hybridized separately per sample (Scenario C)?

4. If metadata column is absent, check the cellranger multi config or the wet-lab protocol documentation.

---

## Orthogonal Barcoding for FDP Analysis

FDP (false discovery proportion) analysis requires an **independent identity label** on each cell that is orthogonal to probe barcode assignment. This lets you define ground truth true positives (DE within probe barcode, across biological conditions) and false positives (DE across probe barcodes, within biological condition).

Valid orthogonal labels:

### CRISPR Guide Assignments
- If cells carry distinct CRISPR guide RNAs and guide identity is captured (e.g., via guide capture in a multiome or separate guide library), you can use guide identity as the orthogonal label
- Each guide acts as an independent cell identity label regardless of which probe barcode the cell is in
- In the paper's VCC dataset, CRISPR guide assignments served this role (Fig. 1G, 1H)

### Hashtag Oligonucleotides (HTOs) / Cell Hashing
- If cells were labeled with distinct HTOs (e.g., from CITE-seq or cell hashing) before probe hybridization, and HTO identity is captured separately, this provides an orthogonal label
- The A375 Flex v2 dataset (Supp. Fig. 6) used HTOs for this purpose

### What does NOT work as orthogonal barcoding
- Cluster identity (clusters are derived from gene expression, which is affected by the artifact)
- Cell cycle phase (same issue — derived from expression)
- Any label that is itself derived from the RNA expression data

---

## Special Case: Superloading / NTC Cells

In the VCC dataset, many cells are **NTC (non-targeting control)** cells — they carry a guide RNA that does not target any gene. When split across all 16 probe barcodes, these cells are biologically identical. Any DE between probe barcodes in NTC cells is **pure technical artifact** with no biological signal.

If your dataset has a similar biologically-identical population split across barcodes (e.g., cells without CRISPR perturbation, untreated controls in all barcodes), use this subset as a negative control:

```r
ntc_cells <- WhichCells(obj, expression = guide_identity == "non-targeting")
obj_ntc <- subset(obj, cells = ntc_cells)
Idents(obj_ntc) <- obj_ntc$probe_barcode

de_ntc <- FindAllMarkers(obj_ntc, only.pos = FALSE, min.pct = 0.1,
                         logfc.threshold = 0, test.use = "wilcox")
# All significant hits here are probe barcode artifacts
```
