#!/usr/bin/env Rscript
# ==============================================================================
# fig2b_stacked_barplot.R
# DE gene counts per probe barcode — A375 FLEXv2 stacked barplot (Figure 2B)
# ==============================================================================
#
# Counts the number of up- and downregulated genes per probe barcode (ligand
# condition) from the one-vs-all FindAllMarkers result for the A375 FLEXv2
# dataset, then plots a stacked barplot.
#
# Thresholds: p_val_adj < 0.05 and |avg_log2FC| >= 0.3.
#
# Each probe barcode corresponds to a distinct ligand treatment condition
# (PBS, ANXA1, GDF15, IFNG, MDK, POSTN, SPP1, TNFA). Unlike the VCC stacked
# barplot (Figure 1E), DE here reflects both biological (ligand treatment) and
# technical (probe barcode) effects.
#
# Y-axis fixed at 0–166 to match Fig 1E (VCC Flex v1) for direct comparison.
#
# INPUT
#   results/DE_A375_flexv2_FindAllMarkers_wilcox_minpct0.1.csv
#     — from 09_compute_a375_de.R
#
# OUTPUT
#   results/figures/fig2b_stacked_barplot_DE_per_probe_barcode_a375.pdf
#
# USAGE
#   conda run -n trekker Rscript scripts/figure2/plot/fig2b_stacked_barplot.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
})

RESULTS_DIR   <- "results"
FIGURES_DIR   <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

INPUT_CSV     <- file.path(RESULTS_DIR, "DE_A375_flexv2_FindAllMarkers_wilcox_minpct0.1.csv")
LOG2FC_CUTOFF <- 0.3
PVAL_CUTOFF   <- 0.05
YLIM_MAX      <- 166   # matches Fig 1E (VCC Flex v1) y-axis maximum

if (!file.exists(INPUT_CSV)) {
  stop("FindAllMarkers CSV not found: ", INPUT_CSV,
       "\nRun 09_compute_a375_de.R first.")
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

bc_order <- c(paste0("A-A", sprintf("%02d", 1:12)),
              paste0("A-B", sprintf("%02d", 1:4)))

counts <- df |>
  group_by(cluster, direction) |>
  summarise(n_genes = n(), .groups = "drop") |>
  complete(cluster, direction, fill = list(n_genes = 0)) |>
  mutate(
    cluster   = factor(cluster,   levels = bc_order),
    direction = factor(direction, levels = c("Downregulated", "Upregulated"))
  )

message("\nDE counts per probe barcode (ligand condition):")
print(as.data.frame(counts))

# ==============================================================================
# Plot
# ==============================================================================

p <- ggplot(counts, aes(x = cluster, y = n_genes, fill = direction)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Upregulated" = "tomato3",
                               "Downregulated" = "steelblue3")) +
  scale_y_continuous(limits = c(0, YLIM_MAX)) +
  labs(x = "Probe barcode (ligand condition)",
       y = "Number of DE genes", fill = NULL) +
  theme_few() +
  theme(
    axis.title        = element_text(size = 24, color = "black"),
    axis.text         = element_text(size = 20, color = "black"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.text       = element_text(size = 18, color = "black")
  )

ggsave(file.path(FIGURES_DIR, "fig2b_stacked_barplot_DE_per_probe_barcode_a375.pdf"),
       p, width = 8.5, height = 5, device = cairo_pdf)
message("Saved: fig2b_stacked_barplot_DE_per_probe_barcode_a375.pdf")

message("Done.")
print(sessionInfo())
