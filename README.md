# flex-technical-variation

Code for reproducing all figures in:

**"Probe barcode-associated technical variation in probe-based single-cell RNA sequencing"**
Jackson A. Weir, Yonit Krebs, Fei Chen — Broad Institute of MIT and Harvard

---

## Overview

This repository contains the scripts to reproduce Figures 1 and 2 from the paper. The analysis characterizes technical gene expression variation associated with probe barcode identity in 10x Genomics Flex (Fixed RNA Profiling) single-cell RNA-seq data, and compares the magnitude of this artifact between Flex v1 and Flex v2.

## Data

All required data files and download instructions are in [`data/README.md`](data/README.md).

**Public data sources used in this paper:**

| Dataset | Source | Figures |
|---------|--------|---------|
| VCC Flex v1 (Arc Institute Virtual Cell Challenge) | [arc-virtual-cell-atlas](https://github.com/ArcInstitute/arc-virtual-cell-atlas) | 1B–1F, 1I, 2D–2E |
| 10x PBMC 320k 16-plex GEM-X Flex | [10x Genomics website](https://www.10xgenomics.com/datasets/320k_Human_PBMCs_Sub_Pool_16-plex_GEM-X_FLEX) | 1I |
| A375 Flex v2 + probe expression caches | [Zenodo 10.5281/zenodo.19363777](https://zenodo.org/records/19363777) | 2B–2F |
| Flex v1 probe set manifest (v1.0.1) | 10x Genomics support site | 2D |
| Flex v2 probe set manifest (v2.0.0) | 10x Genomics support site | 2D |

## Environments

Scripts require either the `trekker` (R 4.3.3) or `mel_spatial` (R 4.5.2) conda environment. Scripts that load `.qs` Seurat objects require `mel_spatial`. Most compute and plot scripts run under `trekker`.

```bash
conda activate trekker   # for most scripts
conda activate mel_spatial  # for scripts loading .qs objects or probe expression caches
```

## Reproducing the figures

### Figure 1 — Flex v1 probe barcode characterization

```
scripts/figure1/
├── 01_setup_vcc_seurat.R                 # annotate VCC Seurat (lane, probe_barcode)
├── 00_prep_violin_data.R                 # extract NTC violin data (mel_spatial)
├── compute/
│   ├── 02_compute_pairwise_de.R          # BC002 vs BC011, lane1 vs lane2 DE
│   ├── 03_compute_findallmarkers_lane1.R # one-vs-all DE per barcode (Lane 1)
│   ├── 04_compute_variance_explained.R   # eta² by probe barcode and lane (NTC)
│   ├── 05_compute_fdp_sweep.R            # FDP sweep BC005 vs BC008
│   ├── 06_setup_pbmc_seurat.R            # build PBMC Seurat
│   └── 07_compute_pbmc_de.R              # BC002 vs BC011 DE in PBMC
└── plot/
    ├── fig1bc_volcanos.R                 # Fig 1B, 1C
    ├── fig1d_violin_top_genes.R          # Fig 1D (requires 00_prep_violin_data.R)
    ├── fig1e_stacked_barplot.R           # Fig 1E
    ├── fig1f_ranked_eta_sq.R             # Fig 1F
    ├── fig1g_volcano_tp_fp.R             # Fig 1G
    ├── fig1h_fdp_sweep.R                 # Fig 1H
    └── fig1i_cross_dataset_scatter.R     # Fig 1I
```

**Run order:**

```bash
# mel_spatial env (loads Seurat .qs objects)
conda run -n mel_spatial Rscript scripts/figure1/01_setup_vcc_seurat.R
conda run -n mel_spatial Rscript scripts/figure1/00_prep_violin_data.R

# trekker env (all other figure 1 scripts)
conda run -n trekker Rscript scripts/figure1/compute/02_compute_pairwise_de.R
conda run -n trekker Rscript scripts/figure1/compute/03_compute_findallmarkers_lane1.R
conda run -n trekker Rscript scripts/figure1/compute/04_compute_variance_explained.R
conda run -n trekker Rscript scripts/figure1/compute/05_compute_fdp_sweep.R
conda run -n trekker Rscript scripts/figure1/compute/06_setup_pbmc_seurat.R
conda run -n trekker Rscript scripts/figure1/compute/07_compute_pbmc_de.R

# Plots (trekker, except fig1d which needs mel_spatial for the prep step)
conda run -n trekker Rscript scripts/figure1/plot/fig1bc_volcanos.R
conda run -n trekker Rscript scripts/figure1/plot/fig1d_violin_top_genes.R
conda run -n trekker Rscript scripts/figure1/plot/fig1e_stacked_barplot.R
conda run -n trekker Rscript scripts/figure1/plot/fig1f_ranked_eta_sq.R
conda run -n trekker Rscript scripts/figure1/plot/fig1g_volcano_tp_fp.R
conda run -n trekker Rscript scripts/figure1/plot/fig1h_fdp_sweep.R
conda run -n trekker Rscript scripts/figure1/plot/fig1i_cross_dataset_scatter.R
```

### Figure 2 — Flex v2 comparison and probe-level mechanism

```
scripts/figure2/
├── 08_setup_a375_seurat.R                       # filter A375 Flex v2 Seurat
├── 00_prep_violin_data_flexv2.R                 # extract A375 violin data (mel_spatial)
├── compute/
│   ├── 09_compute_a375_de.R                     # FindAllMarkers A375 (all barcodes)
│   ├── 10_compute_variance_explained_a375.R     # eta² A375 Flex v2
│   ├── 11_compute_vcc_probe_expression.R        # VCC probe expression cache
│   ├── 12_compute_a375_probe_expression.R       # A375 probe expression cache
│   └── 13_compute_probe_histogram_data.R        # probe counts per gene (Fig 2D inputs)
└── plot/
    ├── fig2b_stacked_barplot.R                  # Fig 2B
    ├── fig2c_violin_top_genes.R                 # Fig 2C (requires 00_prep_violin_data_flexv2.R)
    ├── fig2d_probe_number_histogram.R           # Fig 2D
    ├── fig2e_vcc_probe_expression.R             # Fig 2E
    └── fig2f_a375_probe_expression.R            # Fig 2F
```

**Run order:**

```bash
# mel_spatial env
conda run -n mel_spatial Rscript scripts/figure2/08_setup_a375_seurat.R
conda run -n mel_spatial Rscript scripts/figure2/00_prep_violin_data_flexv2.R
conda run -n mel_spatial Rscript scripts/figure2/compute/11_compute_vcc_probe_expression.R
conda run -n mel_spatial Rscript scripts/figure2/compute/12_compute_a375_probe_expression.R

# trekker env
conda run -n trekker Rscript scripts/figure2/compute/09_compute_a375_de.R
conda run -n trekker Rscript scripts/figure2/compute/10_compute_variance_explained_a375.R
conda run -n trekker Rscript scripts/figure2/compute/13_compute_probe_histogram_data.R

# Plots (trekker, except fig2e/2f which require mel_spatial for Matrix/dittoSeq)
conda run -n trekker Rscript scripts/figure2/plot/fig2b_stacked_barplot.R
conda run -n trekker Rscript scripts/figure2/plot/fig2c_violin_top_genes.R
conda run -n trekker Rscript scripts/figure2/plot/fig2d_probe_number_histogram.R
conda run -n mel_spatial Rscript scripts/figure2/plot/fig2e_vcc_probe_expression.R
conda run -n mel_spatial Rscript scripts/figure2/plot/fig2f_a375_probe_expression.R
```

## Output

All intermediate results are written to `results/`. All figures are saved as Cairo PDF to `results/figures/`.
