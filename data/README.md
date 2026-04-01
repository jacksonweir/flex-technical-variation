# Data

Large data files (Seurat objects, count matrices, expression caches) are not
stored in this repository. This document describes every required file, where
to obtain it, and which figures it feeds.

---

## Data flow overview

```
Public data (download)
    │
    ├─ VCC adata_Training.h5ad ──────────────────────────────► scripts/figure1/01_setup_vcc_seurat.R
    │                                                               │
    │                                                               ▼ results/seurat_vcc_annotated.qs
    │                                                               │
    │                                              ┌────────────────┼──────────────────────────────┐
    │                                              ▼                ▼                              ▼
    │                                    02_compute_pairwise_de  04_compute_variance_explained   00_prep_violin_data.R
    │                                    03_compute_findall...   05_compute_fdp_sweep            (Fig 1D)
    │                                    13_compute_probe_hist.. (Fig 1E)  (Fig 1F,1H,1G)
    │                                    (Fig 2D – VCC side)
    │
    ├─ PBMC filtered_feature_bc_matrix.h5 ──────────────────────► 06_setup_pbmc_seurat.R
    │                                                               │
    │                                                               ▼ 07_compute_pbmc_de.R → Fig 1I
    │
    ├─ seurat_A375_flexv2_filtered.qs (Zenodo) ──────────────────► results/seurat_a375_flexv2_filtered.qs
    │   [place in results/ or run 08_setup_a375_seurat.R]           │
    │                                              ┌─────────────────┼───────────────────┐
    │                                              ▼                 ▼                   ▼
    │                                    09_compute_a375_de   10_compute_var_expl   00_prep_violin_data_flexv2.R
    │                                    (Fig 2B, 2D)          (Fig 2D)              (Fig 2C)
    │
    ├─ vcc_probe_expr_sparse_cache_10k_feat_filt.rds (Zenodo) ──► 11_compute_vcc_probe_expression.R
    │                                                               │
    │                                                               ▼ results/vcc_probe_expression_cache.rds
    │                                                               └──────────────────────────► Fig 2E
    │
    ├─ flex_v2_probe_expr_sparse_cache_10k_feat_filt.rds (Zenodo) ► 12_compute_a375_probe_expression.R
    │                                                               │
    │                                                               ▼ results/a375_probe_expression_cache.rds
    │                                                               └──────────────────────────► Fig 2F
    │
    ├─ vcc_probe_manifest.csv (10x) ──────────────────────────────► 13_compute_probe_histogram_data.R
    └─ a375_v2_probe_manifest.csv (10x) ─────────────────────────► 13_compute_probe_histogram_data.R
                                                                     │
                                                                     ▼ results/probe_counts_vcc.csv
                                                                       results/probe_counts_a375.csv
                                                                     └──────────────────────────► Fig 2D
```

---

## File-by-file details

### Probe manifests

Download and place at the paths below before running scripts.

#### Flex v1 probe set — `data/vcc_probe_manifest.csv`

The VCC dataset was processed with the **Chromium Human Transcriptome Probe Set
v1.0.1** (GRCh38-2020-A). Download from the 10x Genomics support site:

  https://www.10xgenomics.com/support/software/cell-ranger/downloads

Search for "Chromium Human Transcriptome Probe Set v1.0.1". Save the CSV to
`data/vcc_probe_manifest.csv`.

Expected columns: `gene_id, probe_seq, probe_id, included, region`
Expected rows: ~54,580 probes

Used by: `scripts/figure2/compute/13_compute_probe_histogram_data.R` → **Fig 2D**

#### Flex v2 probe set — `data/a375_v2_probe_manifest.csv`

The A375 FLEXv2 dataset was processed with the **Chromium Human Transcriptome
Probe Set v2.0.0** (GRCh38-2024-A). Download from the same 10x Genomics page.

Save to `data/a375_v2_probe_manifest.csv`.

Expected columns: `gene_id, probe_seq, probe_id, included, region, gene_name`

Used by: `scripts/figure2/compute/13_compute_probe_histogram_data.R` → **Fig 2D**

---

### VCC (Virtual Cell Challenge) dataset

**Source:** Arc Institute Virtual Cell Challenge
https://github.com/ArcInstitute/arc-virtual-cell-atlas

Download `adata_Training.h5ad`. Set the path in
`scripts/figure1/01_setup_vcc_seurat.R`:

```r
VCC_H5AD_PATH <- "/path/to/adata_Training.h5ad"
```

This script converts the AnnData object to a Seurat object and adds
lane/probe_barcode annotations, saving `results/seurat_vcc_annotated.qs`.

Used by: all Figure 1 scripts, and Fig 2D (VCC side)

---

### 10x PBMC dataset

**Source:** 10x Genomics
https://www.10xgenomics.com/datasets/320k_Human_PBMCs_Sub_Pool_16-plex_GEM-X_FLEX

Download the `filtered_feature_bc_matrix.h5` file. Set the path in
`scripts/figure1/compute/06_setup_pbmc_seurat.R`:

```r
PBMC_H5_PATH <- "/path/to/filtered_feature_bc_matrix.h5"
```

Used by: `06_setup_pbmc_seurat.R`, `07_compute_pbmc_de.R` → **Fig 1I**

---

### Zenodo deposit — probe barcode technical variation dataset

**DOI:** https://doi.org/10.5281/zenodo.19363777

Three files are available from this deposit. Download them and set paths
in the scripts that use them (see below).

#### `seurat_A375_flexv2_filtered.qs`

Seurat v5 object containing filtered A375 Flex v2 cells, with
`probe_barcode` metadata column (16 barcodes: A-A01 through A-B04).
Gene expression stored as normalized counts.

**Option A** — Place as `results/seurat_a375_flexv2_filtered.qs` to use
directly without running `08_setup_a375_seurat.R`.

**Option B** — Set path in `scripts/figure2/08_setup_a375_seurat.R`:
```r
A375_QS_PATH <- "/path/to/seurat_A375_flexv2_filtered.qs"
```
and run script 08 to apply QC filters and save to `results/`.

Used by: `08_setup_a375_seurat.R` (or placed directly) → scripts 09, 10,
`00_prep_violin_data_flexv2.R` → **Fig 2B, 2C, 2D**

#### `vcc_probe_expr_sparse_cache_10k_feat_filt.rds`

Sparse matrix cache of per-probe (not per-gene) normalized expression for
VCC NTC cells (Flex v1). List with `$data` (dgCMatrix, probes × cells) and
`$barcodes` (probe barcode per cell). Top 10,000 most expressed probes.

Set path in `scripts/figure2/compute/11_compute_vcc_probe_expression.R`:
```r
VCC_PROBE_CACHE <- "/path/to/vcc_probe_expr_sparse_cache_10k_feat_filt.rds"
```

Used by: `11_compute_vcc_probe_expression.R` → `fig2e_vcc_probe_expression.R` → **Fig 2E**

#### `flex_v2_probe_expr_sparse_cache_10k_feat_filt.rds`

Sparse matrix cache of per-probe normalized expression for A375 Flex v2
PBS cells. Same structure as the VCC cache above.

Set path in `scripts/figure2/compute/12_compute_a375_probe_expression.R`:
```r
A375_PROBE_CACHE <- "/path/to/flex_v2_probe_expr_sparse_cache_10k_feat_filt.rds"
```

Used by: `12_compute_a375_probe_expression.R` → `fig2f_a375_probe_expression.R` → **Fig 2F**

---

## Summary: file → figure mapping

| File | Source | Figures |
|------|--------|---------|
| `adata_Training.h5ad` | Arc Institute VCC | 1B, 1C, 1D, 1E, 1F, 1G, 1H, 1I, 2D, 2E |
| `filtered_feature_bc_matrix.h5` | 10x Genomics PBMC | 1I |
| `seurat_A375_flexv2_filtered.qs` | Zenodo | 2B, 2C, 2D |
| `vcc_probe_expr_sparse_cache_10k_feat_filt.rds` | Zenodo | 2E |
| `flex_v2_probe_expr_sparse_cache_10k_feat_filt.rds` | Zenodo | 2F |
| `data/vcc_probe_manifest.csv` | 10x Genomics | 2D |
| `data/a375_v2_probe_manifest.csv` | 10x Genomics | 2D |
