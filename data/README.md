# Data

Large data files are not stored in this repository. Download the files below and set their paths in the relevant scripts.

---

## Zenodo — [10.5281/zenodo.19363777](https://doi.org/10.5281/zenodo.19363777)

| File | Used by | Figures |
|------|---------|---------|
| `seurat_A375_flexv2_filtered.qs` | `scripts/figure2/08_setup_a375_seurat.R` | 2B, 2C, 2D |
| `vcc_probe_expr_sparse_cache_10k_feat_filt.rds` | `scripts/figure2/compute/11_compute_vcc_probe_expression.R` | 2E |
| `flex_v2_probe_expr_sparse_cache_10k_feat_filt.rds` | `scripts/figure2/compute/12_compute_a375_probe_expression.R` | 2F |

## Arc Institute — Virtual Cell Challenge

Download `adata_Training.h5ad` from [arc-virtual-cell-atlas](https://github.com/ArcInstitute/arc-virtual-cell-atlas). Set path in `scripts/figure1/01_setup_vcc_seurat.R`.

Used by: all Figure 1 scripts, Fig 2D, 2E.

## 10x Genomics — PBMC dataset

Download `filtered_feature_bc_matrix.h5` from the [10x Genomics website](https://www.10xgenomics.com/datasets/320k_Human_PBMCs_Sub_Pool_16-plex_GEM-X_FLEX). Set path in `scripts/figure1/compute/06_setup_pbmc_seurat.R`.

Used by: Fig 1I.

## 10x Genomics — Probe set manifests

Download from the [10x Genomics support site](https://www.10xgenomics.com/support/software/cell-ranger/downloads) and save to `data/`:

- **Chromium Human Transcriptome Probe Set v1.1.0 (GRCh38-2024-A)** → `data/vcc_probe_manifest.csv`
- **Chromium Human Transcriptome Probe Set v2.0.0 (GRCh38-2024-A)** → `data/a375_v2_probe_manifest.csv`

Used by: `scripts/figure2/compute/13_compute_probe_histogram_data.R` → Fig 2D.
