# Extracting Probe Barcode Identity from Cell Barcodes (Flex v1)

## Background

In 10x Genomics Flex v1, the probe barcode identity is encoded directly in the cell barcode.
Each cell barcode is 16 bp; the **last 8 bp of that 16 bp sequence identify which of the 16
probe barcodes (BC001–BC016) that cell came from.**

When cellranger multi processes a Flex v1 library, it assigns cells to per-sample outputs
based on this embedded sequence. If you loaded data directly from the per-sample outputs and
merged the objects yourself, the probe barcode identity may already be present as a metadata
column. If you loaded a combined matrix (e.g., from cellranger aggr or a merged Seurat object)
without per-sample metadata, use the code below to recover probe barcode assignments from the
cell barcode sequences.

This approach is implemented in `scripts/figure1/01_setup_vcc_seurat.R` in this repository.

---

## Probe Barcode Sequence Whitelist (Flex v1)

The 8 bp sequences embedded in cell barcodes for all 16 Flex v1 probe barcodes:

| Probe barcode | 8 bp sequence |
|--------------|---------------|
| BC001 | ACTTTAGG |
| BC002 | AACGGGAA |
| BC003 | AGTAGGCT |
| BC004 | ATGTTGAC |
| BC005 | ACAGACCT |
| BC006 | ATCCCAAC |
| BC007 | AAGTAGAG |
| BC008 | AGCTGTGA |
| BC009 | ACAGTCTG |
| BC010 | AGTGAGTG |
| BC011 | AGAGGCAA |
| BC012 | ACTACTCA |
| BC013 | ATACGTCA |
| BC014 | ATCATGTG |
| BC015 | AACGCCGA |
| BC016 | ATTCGGTT |

Source: 10x Genomics `probe-barcodes-fixed-rna-profiling.txt`.

---

## R Code

```r
library(Seurat)
library(tibble)

# Probe barcode whitelist — Flex v1, BC001–BC016
probe_barcode_ref <- tibble(
  probe_barcode = paste0("BC", sprintf("%03d", 1:16)),
  sequence      = c("ACTTTAGG", "AACGGGAA", "AGTAGGCT", "ATGTTGAC",
                    "ACAGACCT", "ATCCCAAC", "AAGTAGAG", "AGCTGTGA",
                    "ACAGTCTG", "AGTGAGTG", "AGAGGCAA", "ACTACTCA",
                    "ATACGTCA", "ATCATGTG", "AACGCCGA", "ATTCGGTT")
)

# Cell barcodes follow the format: <16 bp sequence>-<numeric suffix>
# e.g. "ACGTACGTACGTACGT-1"
# Strip the suffix to get the 16 bp portion, then take the last 8 characters.
bc_16bp <- sub("-.*$", "", colnames(seurat_obj))
seurat_obj$probe_barcode_seq <- substr(bc_16bp, nchar(bc_16bp) - 7, nchar(bc_16bp))

# Map the 8 bp sequence to a probe barcode ID (BC001–BC016)
seurat_obj$probe_barcode <- probe_barcode_ref$probe_barcode[
  match(seurat_obj$probe_barcode_seq, probe_barcode_ref$sequence)
]

# Check how many cells were matched
n_matched   <- sum(!is.na(seurat_obj$probe_barcode))
n_unmatched <- sum(is.na(seurat_obj$probe_barcode))
message(sprintf("Matched: %d cells  |  Unmatched: %d cells", n_matched, n_unmatched))

# Distribution across barcodes
print(table(seurat_obj$probe_barcode))
```

---

## Expected Output

All (or nearly all) cells should match one of the 16 probe barcodes. A small number of
unmatched cells may occur due to cells near quality thresholds or edge cases in barcode
correction. If a large fraction of cells are unmatched, check:

1. **Is this actually a Flex v1 run?** Flex v2 cell barcodes use a different format and
   the last 8 bp will not match the Flex v1 whitelist.
2. **Are barcodes in the expected format?** Run `head(colnames(seurat_obj))` to confirm
   they follow `<16bp>-<suffix>`. Some preprocessing pipelines truncate barcodes or drop
   the suffix entirely.
3. **Were multiple lanes merged?** Cells from different lanes may carry different numeric
   suffixes (-1, -2, -3). The extraction still works because it operates on the 16 bp
   portion, not the suffix.

---

## Note: Flex v2

For Flex v2, the plate-well identity (e.g., `A-A01`) is assigned during the separate oligo
barcoding step and is provided by cellranger multi directly — it does not need to be extracted
from cell barcodes. The `probe_barcodes_v2.txt` file in this repository contains the full
Flex v2 barcode whitelist with plate-well annotations. Consult your cellranger multi config
or the per-sample output directories to confirm which plate-well corresponds to which biological
condition.
