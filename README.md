# flex-technical-variation

Code for reproducing all figures in:

**"Sample barcoding-associated technical variation in probe-based single-cell RNA sequencing"**

In 10x Genomics Flex v1, probe set barcode identity drives substantial technical variation in gene expression — reproducible across lanes and datasets — that generates a high rate of false positive DE genes when barcodes are confounded with biological sample identity. The Flex v2 assay, which decouples sample barcoding from probe hybridization, significantly reduces this artifact.

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
| Flex v1 probe set manifest (v1.1.0) | 10x Genomics support site | 2D |
| Flex v2 probe set manifest (v2.0.0) | 10x Genomics support site | 2D |

## Reproducing the figures

### Figure 1 — Flex v1 probe barcode characterization

```
scripts/figure1/
├── 01_setup_vcc_seurat.R                 # annotate VCC Seurat (lane, probe_barcode)
├── 00_prep_violin_data.R                 # extract NTC violin data
├── compute/
│   ├── 02_compute_pairwise_de.R          # BC002 vs BC011, lane1 vs lane2 DE
│   ├── 03_compute_findallmarkers_lane1.R # one-vs-all DE per barcode (Lane 1)
│   ├── 04_compute_variance_explained.R   # eta² by probe barcode and lane (NTC)
│   ├── 05_compute_fdp_sweep.R            # FDP sweep BC005 vs BC008
│   ├── 06_setup_pbmc_seurat.R            # build PBMC Seurat
│   └── 07_compute_pbmc_de.R              # BC002 vs BC011 DE in PBMC
└── plot/
    ├── fig1bc_volcanos.R                 # Fig 1B, 1C
    ├── fig1d_violin_top_genes.R          # Fig 1D
    ├── fig1e_stacked_barplot.R           # Fig 1E
    ├── fig1f_ranked_eta_sq.R             # Fig 1F
    ├── fig1g_volcano_tp_fp.R             # Fig 1G
    ├── fig1h_fdp_sweep.R                 # Fig 1H
    └── fig1i_cross_dataset_scatter.R     # Fig 1I
```

### Figure 2 — Flex v2 comparison and probe-level mechanism

```
scripts/figure2/
├── 08_setup_a375_seurat.R                       # filter A375 Flex v2 Seurat
├── 00_prep_violin_data_flexv2.R                 # extract A375 violin data
├── compute/
│   ├── 09_compute_a375_de.R                     # FindAllMarkers A375 (all barcodes)
│   ├── 10_compute_variance_explained_a375.R     # eta² A375 Flex v2
│   ├── 11_compute_vcc_probe_expression.R        # VCC probe expression cache
│   ├── 12_compute_a375_probe_expression.R       # A375 probe expression cache
│   └── 13_compute_probe_histogram_data.R        # probe counts per gene (Fig 2D inputs)
└── plot/
    ├── fig2b_stacked_barplot.R                  # Fig 2B
    ├── fig2c_violin_top_genes.R                 # Fig 2C
    ├── fig2d_probe_number_histogram.R           # Fig 2D
    ├── fig2e_vcc_probe_expression.R             # Fig 2E
    └── fig2f_a375_probe_expression.R            # Fig 2F
```

## Agentic Skill — `flex-artifact-skill`

The `flex-artifact-skill/` directory contains a structured skill for AI coding agents (Claude Code, GitHub Copilot, Codex, etc.). The skill provides agents with the context they need to:

- Understand what probe barcode-associated technical variation is and how it was characterized in this paper
- Apply the same analytical approaches (η², FDP sweep, FindAllMarkers) to a user's own Flex dataset
- Assess whether a user's differential expression results are likely inflated by probe barcode artifacts, given their experimental design
- Interpret findings in light of the paper's benchmarks (FDP curves, sentinel genes, cross-dataset reproducibility)
- Recommend appropriate mitigation strategies

The skill follows the [Claude skill format](https://github.com/K-Dense-AI/claude-scientific-skills): a concise `SKILL.md` that an agent loads when the topic is relevant, with detailed reference files in `references/` that are read on demand to avoid bloating the context window.

```
flex-artifact-skill/
├── SKILL.md                          # Overview, 7-step workflow, quick reference tables
└── references/
    ├── experimental-design.md        # Flex v1 vs v2 chemistry, scenarios, how to identify your design
    ├── analysis-code.md              # Full R/Python code: η², FDP, DE, all visualizations
    └── interpretation.md             # Risk table, paper benchmarks, sentinel genes, mitigation strategies
```

The supplementary DE tables (`supplementary_tables/`) are also included for direct cross-referencing: overlap between a user's DE genes and these reference tables is strong evidence of probe barcode confounding.
