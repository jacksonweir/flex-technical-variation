#!/usr/bin/env Rscript
# ==============================================================================
# fig1e_stacked_barplot.R
# DE gene counts per probe barcode — stacked barplot (Figure 1E)
# ==============================================================================
#
# Loads the FindAllMarkers result for Lane 1 and counts the number of
# upregulated and downregulated genes per probe barcode after applying
# p_val_adj < 0.05 and |log2FC| >= 0.3 thresholds.
#
# Each probe barcode was compared against all other probe barcodes pooled
# (one-vs-all). Bars are ordered BC001–BC016.
#
# INPUT
#   results/DE_FindAllMarkers_lane1_wilcox.csv  — from 03_compute_findallmarkers_lane1.R
#
# OUTPUT
#   results/figures/fig1e_stacked_barplot_DE_per_probe_barcode.pdf
#
# USAGE
#   Rscript scripts/figure1/plot/fig1e_stacked_barplot.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

LOG2FC_CUTOFF <- 0.3
PVAL_CUTOFF   <- 0.05

INPUT_CSV <- file.path(RESULTS_DIR, "DE_FindAllMarkers_lane1_wilcox.csv")

if (!file.exists(INPUT_CSV)) {
  stop("FindAllMarkers CSV not found: ", INPUT_CSV,
       "\nRun 03_compute_findallmarkers_lane1.R first.")
}

# ==============================================================================
# Load, filter, and count per barcode
# ==============================================================================

df <- read.csv(INPUT_CSV) |>
  filter(!is.na(avg_log2FC), !is.na(p_val_adj)) |>
  filter(abs(avg_log2FC) >= LOG2FC_CUTOFF, p_val_adj < PVAL_CUTOFF) |>
  mutate(direction = case_when(
    avg_log2FC >= LOG2FC_CUTOFF  ~ "Upregulated",
    avg_log2FC <= -LOG2FC_CUTOFF ~ "Downregulated"
  )) |>
  filter(!is.na(direction))

message(sprintf("Significant genes: %d  (Up: %d  Down: %d)",
  nrow(df), sum(df$direction == "Upregulated"), sum(df$direction == "Downregulated")))

counts <- df |>
  group_by(cluster, direction) |>
  summarise(n_genes = n(), .groups = "drop") |>
  complete(cluster, direction, fill = list(n_genes = 0)) |>
  mutate(
    cluster   = factor(cluster, levels = paste0("BC", sprintf("%03d", 1:16))),
    direction = factor(direction, levels = c("Downregulated", "Upregulated"))
  )

# ==============================================================================
# Plot
# ==============================================================================

p <- ggplot(counts, aes(x = cluster, y = n_genes, fill = direction)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Upregulated" = "tomato3",
                               "Downregulated" = "steelblue3")) +
  labs(x = "Probe barcode", y = "Number of DE genes", fill = NULL) +
  theme_few() +
  theme(
    axis.title        = element_text(size = 20, color = "black"),
    axis.text         = element_text(size = 16, color = "black"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.text       = element_text(size = 15, color = "black")
  )

ggsave(file.path(FIGURES_DIR, "fig1e_stacked_barplot_DE_per_probe_barcode.pdf"),
       p, width = 8, height = 5, device = cairo_pdf)
message("Saved: fig1e_stacked_barplot_DE_per_probe_barcode.pdf")

message("Done.")
print(sessionInfo())
