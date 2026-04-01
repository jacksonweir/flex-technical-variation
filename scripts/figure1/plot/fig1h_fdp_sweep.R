#!/usr/bin/env Rscript
# ==============================================================================
# fig1h_fdp_sweep.R
# False discovery proportion vs |log2FC| threshold sweep (Figure 1H)
# ==============================================================================
#
# Plots FDP = FP / (FP + TP) as a function of |log2FC| threshold for the
# BC005 vs BC008 barcode pair (representing the average technical effect among
# the 16 probe barcodes). Each grey line is one target gene; the black line
# and points show the median across all genes.
#
# INPUT
#   results/fdp_sweep_BC005_vs_BC008/metrics_all_genes.csv
#     — from 05_compute_fdp_sweep.R
#
# OUTPUT
#   results/figures/fig1h_fdp_sweep_BC005_vs_BC008.pdf
#
# USAGE
#   Rscript scripts/figure1/plot/fig1h_fdp_sweep.R
#
# PACKAGES
#   ggplot2, dplyr
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

INPUT_CSV <- file.path(RESULTS_DIR, "fdp_sweep_BC005_vs_BC008", "metrics_all_genes.csv")

if (!file.exists(INPUT_CSV)) {
  stop("FDP sweep metrics not found: ", INPUT_CSV,
       "\nRun 05_compute_fdp_sweep.R first.")
}

# ==============================================================================
# Load and compute median FDP at each threshold
# ==============================================================================

df_all <- read.csv(INPUT_CSV)

df_median <- df_all |>
  group_by(log2fc_cutoff) |>
  summarise(med_fdp = median(FDP, na.rm = TRUE), .groups = "drop")

message(sprintf("Genes: %d  |  Thresholds: %d",
  length(unique(df_all$target_gene)),
  length(unique(df_all$log2fc_cutoff))))
message(sprintf("Median FDP at |log2FC| >= 0.3: %.2f",
  df_median$med_fdp[df_median$log2fc_cutoff == 0.3]))

# ==============================================================================
# Plot
# ==============================================================================

p <- ggplot() +
  geom_line(
    data      = df_all,
    aes(x = log2fc_cutoff, y = FDP, group = target_gene),
    linewidth = 0.35, color = "grey50", alpha = 0.10
  ) +
  geom_line(
    data      = df_median,
    aes(x = log2fc_cutoff, y = med_fdp),
    linewidth = 1.6, color = "black"
  ) +
  geom_point(
    data  = df_median,
    aes(x = log2fc_cutoff, y = med_fdp),
    size  = 2.2, color = "black", alpha = 0.95
  ) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "|log2FC| cutoff", y = "False Discovery Proportion") +
  theme_classic() +
  theme(
    legend.position   = "none",
    text              = element_text(color = "black"),
    axis.title        = element_text(size = 20, color = "black"),
    axis.text         = element_text(size = 16, color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "fig1h_fdp_sweep_BC005_vs_BC008.pdf"),
       p, width = 7, height = 5, device = cairo_pdf)
message("Saved: fig1h_fdp_sweep_BC005_vs_BC008.pdf")

message("Done.")
print(sessionInfo())
