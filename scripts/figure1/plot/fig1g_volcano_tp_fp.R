#!/usr/bin/env Rscript
# ==============================================================================
# fig1g_volcano_tp_fp.R
# TP / FP volcano for INSIG1 — BC005 vs BC008 confounding (Figure 1G)
# ==============================================================================
#
# Loads the ground truth DE and confounded DE for the gene INSIG1 from the
# BC005 vs BC008 FDP sweep. Classifies each gene in the confounded analysis as:
#
#   True Positive  — significant in both ground truth and confounded analyses
#   False Positive — significant in confounded but not in ground truth
#   Not Significant — not significant in confounded analysis
#
# Illustrates that false positive genes are indistinguishable from true
# positives by effect size or p-value.
#
# Significance threshold: p_val_adj < 0.05 and |log2FC| >= 0.3.
#
# INPUT
#   results/fdp_sweep_BC005_vs_BC008/de_truth_INSIG1.csv
#   results/fdp_sweep_BC005_vs_BC008/de_confounded_INSIG1.csv
#     — both from 05_compute_fdp_sweep.R
#
# OUTPUT
#   results/figures/fig1g_volcano_INSIG1_TP_FP.pdf
#
# USAGE
#   Rscript scripts/figure1/plot/fig1g_volcano_tp_fp.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
SWEEP_DIR   <- file.path(RESULTS_DIR, "fdp_sweep_BC005_vs_BC008")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

PVAL_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 0.3
MAX_LABELS    <- 10
TARGET_GENE   <- "INSIG1"

TRUTH_CSV <- file.path(SWEEP_DIR, paste0("de_truth_",      TARGET_GENE, ".csv"))
CONF_CSV  <- file.path(SWEEP_DIR, paste0("de_confounded_", TARGET_GENE, ".csv"))

for (f in c(TRUTH_CSV, CONF_CSV)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f, "\nRun 05_compute_fdp_sweep.R first.")
  }
}

# ==============================================================================
# Load and classify
# ==============================================================================

df_truth <- read.csv(TRUTH_CSV)
df_conf  <- read.csv(CONF_CSV)

truth_de_genes <- df_truth$gene[
  df_truth$p_val_adj < PVAL_CUTOFF & abs(df_truth$avg_log2FC) >= LOG2FC_CUTOFF
]
message(sprintf("Ground truth DE genes for %s: %d", TARGET_GENE, length(truth_de_genes)))

df_plot <- df_conf %>%
  mutate(
    p_val_adj  = ifelse(is.na(p_val_adj), 1, p_val_adj),
    log10_pval = -log10(p_val_adj + 1e-300),
    is_sig     = p_val_adj < PVAL_CUTOFF & abs(avg_log2FC) >= LOG2FC_CUTOFF,
    reg_dir    = case_when(
      !is_sig                              ~ "Not Significant",
      is_sig & gene %in% truth_de_genes   ~ "True Positive",
      is_sig & !gene %in% truth_de_genes  ~ "False Positive"
    )
  )

n_tp <- sum(df_plot$reg_dir == "True Positive")
n_fp <- sum(df_plot$reg_dir == "False Positive")
message(sprintf("Confounded analysis — TP: %d  FP: %d  FDP: %.1f%%",
                n_tp, n_fp, 100 * n_fp / max(n_tp + n_fp, 1)))

# Label top significant genes by -log10(p) (split evenly TP/FP)
df_label <- df_plot %>%
  filter(reg_dir %in% c("True Positive", "False Positive")) %>%
  arrange(desc(log10_pval)) %>%
  group_by(reg_dir) %>%
  slice_head(n = ceiling(MAX_LABELS / 2)) %>%
  ungroup()

df_plot$label  <- ifelse(df_plot$gene %in% df_label$gene, df_plot$gene, NA_character_)
df_plot$reg_dir <- factor(df_plot$reg_dir,
  levels = c("True Positive", "False Positive", "Not Significant"))

# ==============================================================================
# Plot
# ==============================================================================

p <- ggplot(df_plot, aes(x = avg_log2FC, y = log10_pval)) +
  geom_point(aes(color = reg_dir, shape = reg_dir), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c(
    "True Positive"   = "#3CB371",
    "False Positive"  = "tomato3",
    "Not Significant" = "grey"
  )) +
  scale_shape_manual(values = c(
    "True Positive"   = 16,
    "False Positive"  = 8,
    "Not Significant" = 16
  )) +
  geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(PVAL_CUTOFF),
             linetype = "dashed", color = "black") +
  geom_text_repel(
    aes(label = label),
    size = 3, force = 0.5, max.overlaps = MAX_LABELS,
    box.padding = 0.1, point.padding = 0.1, segment.size = 0.1,
    na.rm = TRUE, color = "black", segment.color = "black"
  ) +
  labs(
    x = sprintf("Average log2 Fold Change (%s_BC008 / NTC_BC005)", TARGET_GENE),
    y = "-log10(adjusted p-value)"
  ) +
  theme_classic() +
  theme(
    text              = element_text(color = "black"),
    axis.title        = element_text(size = 17, color = "black"),
    axis.text         = element_text(size = 14, color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "none"
  )

ggsave(file.path(FIGURES_DIR, "fig1g_volcano_INSIG1_TP_FP.pdf"),
       p, width = 7, height = 5, device = cairo_pdf)
message("Saved: fig1g_volcano_INSIG1_TP_FP.pdf")

message("Done.")
print(sessionInfo())
