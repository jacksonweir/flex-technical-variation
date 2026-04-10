# VCC Reference Comparison and Confounded Design Check

## Purpose

This reference describes two complementary analyses for checking whether probe barcode
technical artifacts are present in your DE results:

1. **One-vs-all scatter plots per barcode (primary)** — Run FindAllMarkers on your data
   with probe barcode as the identity. For each barcode, plot your avg_log2FC against the
   VCC reference avg_log2FC for the same barcode. Since both use the same one-vs-all
   approach, the fold-changes are directly comparable. A positive Spearman ρ across barcodes
   is strong evidence that the same artifact is present in your data.

2. **Pairwise comparison for specific barcode pairs (advanced)** — If your experimental
   design is confounded (one biological condition on one barcode), run pairwise DE for that
   specific barcode pair and compare to the VCC reference for the same pair. Useful for
   quantifying artifact overlap in a specific comparison of interest.

---

## Reference Files

Tables 1 and 2 are available in `references/` in this skill directory:

- `references/Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv` — **VCC Flex v1 dataset,
  all cells, Lane 1, 16 probe barcodes.** This is the primary reference. VCC NTC cells are
  biologically identical across all barcodes, so any DE between barcodes is unambiguously
  technical artifact.

- `references/Supplementary_Table_2_PBMC_Flexv1_FindAllMarkers.csv` — **10x PBMC 320k
  dataset, all cells, 16 probe barcodes.** No CRISPR perturbations; DE between barcodes is
  unambiguously technical. Use as a second independent confirmation.

Both tables have columns: `p_val`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val_adj`,
`cluster` (probe barcode ID, e.g., `"BC002"`), `gene`. A positive `avg_log2FC` in the BC002
cluster means that gene is elevated in BC002 relative to all other barcodes combined.

---

## Step 1 — Run FindAllMarkers on Your Data

Set probe barcode as the active identity and run FindAllMarkers one-vs-all. This mirrors
exactly how the VCC reference (Table 1) was generated, making the fold-changes directly
comparable.

```r
library(Seurat)

Idents(seurat_obj) <- seurat_obj$probe_barcode

user_fam <- FindAllMarkers(
  seurat_obj,
  only.pos        = FALSE,
  min.pct         = 0.1,
  logfc.threshold = 0,
  test.use        = "wilcox",
  verbose         = TRUE
)
user_fam$p_val_adj <- ifelse(is.na(user_fam$p_val_adj), 1, user_fam$p_val_adj)
```

---

## Step 2 — Load the VCC Reference

```r
library(dplyr)

vcc_table1 <- read.csv(
  "references/Supplementary_Table_1_VCC_Flexv1_FindAllMarkers.csv",
  stringsAsFactors = FALSE
)

PVAL_CUTOFF     <- 0.05
FC_CUTOFF       <- 0.3
RECURRENT_GENES <- c("SRP9", "SPCS1", "ISCU", "RRM2", "SLC2A3", "ITSN1", "PHF6", "ADCY2")
```

---

## Step 3 — Per-Barcode One-vs-All Scatter Plots

For each probe barcode, merge your FindAllMarkers results with the VCC reference for that
same barcode cluster and plot avg_log2FC against each other. Because both use the same
one-vs-all framework, the fold-changes are on the same scale.

```r
library(ggplot2)
library(ggrepel)

# Barcodes present in both your data and the VCC reference
barcodes_to_plot <- intersect(
  unique(user_fam$cluster),
  unique(vcc_table1$cluster)
)

plots <- lapply(barcodes_to_plot, function(bc) {

  df_user <- user_fam |>
    filter(cluster == bc) |>
    select(gene, avg_log2FC_user = avg_log2FC, p_val_adj_user = p_val_adj)

  df_vcc <- vcc_table1 |>
    filter(cluster == bc) |>
    select(gene, avg_log2FC_vcc = avg_log2FC, p_val_adj_vcc = p_val_adj)

  df <- merge(df_user, df_vcc, by = "gene")
  if (nrow(df) < 10) return(NULL)

  df$sig_user <- df$p_val_adj_user < PVAL_CUTOFF & abs(df$avg_log2FC_user) >= FC_CUTOFF
  df$sig_vcc  <- df$p_val_adj_vcc  < PVAL_CUTOFF & abs(df$avg_log2FC_vcc)  >= FC_CUTOFF

  df$category <- with(df, case_when(
    sig_user & sig_vcc ~ "Significant in both",
    sig_user           ~ "User only",
    sig_vcc            ~ "VCC only",
    TRUE               ~ "Not significant"
  ))

  # Label genes significant in both OR recurrently affected
  df$label <- ifelse(
    df$category == "Significant in both" | df$gene %in% RECURRENT_GENES,
    df$gene, NA_character_
  )

  cor_r <- cor.test(df$avg_log2FC_user, df$avg_log2FC_vcc, method = "spearman")$estimate
  n_both <- sum(df$category == "Significant in both")

  ggplot(df, aes(x = avg_log2FC_vcc, y = avg_log2FC_user, color = category)) +
    geom_point(alpha = 0.7, size = 1.4, shape = 16) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey60", linewidth = 0.4) +
    scale_color_manual(values = c(
      "Significant in both" = "tomato3",
      "User only"           = "steelblue3",
      "VCC only"            = "orange3",
      "Not significant"     = "grey80"
    )) +
    geom_text_repel(
      aes(label = label),
      size = 2.5, max.overlaps = 8, box.padding = 0.3, point.padding = 0.1,
      color = "black", segment.color = "black", segment.size = 0.2, na.rm = TRUE
    ) +
    labs(
      x     = sprintf("VCC avg log2FC\n(%s vs all)", bc),
      y     = sprintf("Your data avg log2FC\n(%s vs all)", bc),
      title = sprintf("%s  |  ρ = %.2f  |  both sig: %d", bc, cor_r, n_both),
      color = NULL
    ) +
    theme_classic() +
    theme(
      axis.title       = element_text(size = 10, color = "black"),
      axis.text        = element_text(size = 9,  color = "black"),
      plot.title       = element_text(size = 9,  color = "black"),
      axis.ticks       = element_line(color = "black"),
      axis.line        = element_line(color = "black", linewidth = 0.5),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
})

plots <- Filter(Negate(is.null), plots)

# Display as a grid (requires patchwork)
library(patchwork)
wrap_plots(plots, ncol = 4)
```

### Summary statistics across all barcodes

```r
# Spearman ρ per barcode and count of genes significant in both
summary_df <- bind_rows(lapply(barcodes_to_plot, function(bc) {
  df_user <- user_fam |> filter(cluster == bc) |>
    select(gene, avg_log2FC_user = avg_log2FC, p_val_adj_user = p_val_adj)
  df_vcc  <- vcc_table1 |> filter(cluster == bc) |>
    select(gene, avg_log2FC_vcc = avg_log2FC, p_val_adj_vcc = p_val_adj)
  df <- merge(df_user, df_vcc, by = "gene")
  if (nrow(df) < 10) return(NULL)
  cor_r <- cor.test(df$avg_log2FC_user, df$avg_log2FC_vcc, method = "spearman")$estimate
  n_both <- sum(
    df$p_val_adj_user < PVAL_CUTOFF & abs(df$avg_log2FC_user) >= FC_CUTOFF &
    df$p_val_adj_vcc  < PVAL_CUTOFF & abs(df$avg_log2FC_vcc)  >= FC_CUTOFF,
    na.rm = TRUE
  )
  data.frame(barcode = bc, spearman_r = round(cor_r, 3), n_sig_both = n_both)
}))
print(summary_df)
```

---

## Step 4 (Advanced) — Pairwise Comparison for a Specific Barcode Pair

Use this when you want a direct pairwise comparison for a specific pair of barcodes (e.g.,
the two barcodes that correspond to your treatment vs control in a confounded design).
Because Table 1 contains one-vs-all results, this step approximates a pairwise VCC log2FC
by treating genes elevated in one barcode as positive and genes elevated in the other as
negative.

```r
library(Seurat)

BARCODE_HIGH <- "BC002"  # The barcode for condition A (edit to match your design)
BARCODE_LOW  <- "BC011"  # The barcode for condition B

Idents(seurat_obj) <- seurat_obj$probe_barcode

user_de <- FindMarkers(
  seurat_obj,
  ident.1         = BARCODE_HIGH,
  ident.2         = BARCODE_LOW,
  test.use        = "wilcox",
  min.pct         = 0.1,
  logfc.threshold = 0,
  verbose         = FALSE
)
user_de$gene      <- rownames(user_de)
user_de$p_val_adj <- ifelse(is.na(user_de$p_val_adj), 1, user_de$p_val_adj)

# If your design is confounded (probe barcode = sample ID), pass your existing
# sample-vs-sample DE results here instead of re-running FindMarkers.
```

### Reconstruct approximate pairwise VCC log2FC from Table 1

```r
# Genes elevated in BARCODE_HIGH → positive direction
vcc_high <- vcc_table1 |>
  filter(cluster == BARCODE_HIGH, p_val_adj < PVAL_CUTOFF, avg_log2FC > 0) |>
  select(gene, avg_log2FC_vcc = avg_log2FC, p_val_adj_vcc = p_val_adj)

# Genes elevated in BARCODE_LOW → appear DOWN in BARCODE_HIGH, so negate
vcc_low <- vcc_table1 |>
  filter(cluster == BARCODE_LOW, p_val_adj < PVAL_CUTOFF, avg_log2FC > 0) |>
  select(gene, avg_log2FC_vcc = avg_log2FC, p_val_adj_vcc = p_val_adj) |>
  mutate(avg_log2FC_vcc = -avg_log2FC_vcc)

# Combine; prioritise BARCODE_HIGH if a gene appears in both clusters
vcc_ref <- bind_rows(vcc_high, vcc_low) |>
  group_by(gene) |> slice(1) |> ungroup()
```

### Merge, classify, and summarise

```r
df <- merge(
  user_de[, c("gene", "avg_log2FC", "p_val_adj")],
  vcc_ref,
  by = "gene", all.x = TRUE
)
df$p_val_adj[is.na(df$p_val_adj)] <- 1

df$category <- with(df, case_when(
  p_val_adj < PVAL_CUTOFF & !is.na(p_val_adj_vcc) & p_val_adj_vcc < PVAL_CUTOFF
    ~ "Significant in both (likely artifact)",
  p_val_adj < PVAL_CUTOFF & (is.na(p_val_adj_vcc) | p_val_adj_vcc >= PVAL_CUTOFF)
    ~ "User DE only",
  !is.na(p_val_adj_vcc) & p_val_adj_vcc < PVAL_CUTOFF & p_val_adj >= PVAL_CUTOFF
    ~ "VCC artifact only",
  TRUE ~ "Not significant"
))

your_sig_genes <- df$gene[df$p_val_adj < PVAL_CUTOFF & abs(df$avg_log2FC) >= FC_CUTOFF]

cat("=== Confounded Design Check ===\n")
cat(sprintf("Your DE genes at |log2FC| >= %.1f: %d\n", FC_CUTOFF, length(your_sig_genes)))
cat(sprintf("Of these, significant in VCC reference: %d\n",
    sum(df$category == "Significant in both (likely artifact)" &
        abs(df$avg_log2FC) >= FC_CUTOFF, na.rm = TRUE)))
cat(sprintf("Recurrently affected genes in your DE: %s\n",
    paste(intersect(RECURRENT_GENES, your_sig_genes), collapse = ", ") |>
    (\(x) if (nchar(x) == 0) "none" else x)()))
```

### Scatter plot

```r
df_plot <- df[!is.na(df$avg_log2FC_vcc), ]

ggplot(df_plot, aes(x = avg_log2FC, y = avg_log2FC_vcc, color = category)) +
  geom_point(alpha = 0.75, size = 1.5, shape = 16) +
  geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF),
             linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0,
             linetype = "solid", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c(
    "Significant in both (likely artifact)" = "tomato3",
    "User DE only"                          = "steelblue3",
    "VCC artifact only"                     = "orange3",
    "Not significant"                       = "grey80"
  ), drop = FALSE) +
  geom_text_repel(
    data = df_plot[df_plot$category == "Significant in both (likely artifact)", ],
    aes(label = gene),
    size = 3, max.overlaps = 10, box.padding = 0.3, point.padding = 0.1,
    color = "black", segment.color = "black", segment.size = 0.2, na.rm = TRUE
  ) +
  labs(
    x       = sprintf("Your data: avg log2FC (%s vs %s)", BARCODE_HIGH, BARCODE_LOW),
    y       = sprintf("VCC reference: approx log2FC (%s vs %s)", BARCODE_HIGH, BARCODE_LOW),
    color   = NULL,
    caption = "VCC log2FC approximated from one-vs-all FindAllMarkers (Table 1)"
  ) +
  theme_classic() +
  theme(
    axis.title        = element_text(size = 14, color = "black"),
    axis.text         = element_text(size = 12, color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.text       = element_text(size = 10, color = "black")
  )

cor_result <- cor.test(df_plot$avg_log2FC, df_plot$avg_log2FC_vcc, method = "spearman")
cat(sprintf("Spearman r = %.3f  (p = %.2e)\n", cor_result$estimate, cor_result$p.value))
```

---

## Interpreting the Results

**One-vs-all scatter plots (Step 3):**
- A positive Spearman ρ for a given barcode (particularly ρ > 0.3) means your data shows
  the same directional gene shifts for that barcode as in the VCC reference — evidence the
  probe barcode artifact is present.
- Consistent positive ρ across most or all 16 barcodes is strong overall evidence of
  artifact. If only one or two barcodes show correlation, the effect may be concentrated
  in specific barcode pairs.
- Genes labeled red (significant in both) are your highest-priority artifact candidates for
  that barcode.

**Confounded design:**
- If one biological condition is entirely on one probe barcode, probe barcode artifact is
  inseparable from biological DE in your results.
- The recurrently affected genes (SRP9, SPCS1, ISCU, RRM2, SLC2A3, ITSN1, PHF6, ADCY2)
  appearing in your DE list are almost certainly technical artifacts.
- A high Spearman ρ (> 0.4) in the pairwise plot (Step 4) and many genes significant in
  both datasets indicates heavy artifact contamination. Consult `references/interpretation.md`
  for mitigation strategies.

**Using Table 2 (PBMC) as a second check:**
- Repeat Step 3 using `Supplementary_Table_2_PBMC_Flexv1_FindAllMarkers.csv` instead of
  Table 1. Genes consistently flagged by both VCC and PBMC references are highly reliable
  artifact indicators.
