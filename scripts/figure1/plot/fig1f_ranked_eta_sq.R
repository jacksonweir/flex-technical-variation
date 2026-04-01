#!/usr/bin/env Rscript
# ==============================================================================
# fig1f_ranked_eta_sq.R
# Ranked gene plot by probe barcode eta-squared (Figure 1F)
# ==============================================================================
#
# Genes are ranked by eta-squared (variance explained by probe barcode identity)
# and plotted on a scatter. Genes with eta² >= 1% are colored tomato3; the rest
# grey. The top 20 genes above the 1% threshold are labelled.
#
# INPUT
#   results/variance_explained_probe_barcode_NTC.csv  — from 04_compute_variance_explained.R
#
# OUTPUT
#   results/figures/fig1f_ranked_eta_sq_probe_barcode.pdf
#
# USAGE
#   Rscript scripts/figure1/plot/fig1f_ranked_eta_sq.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  library(dplyr)
  library(scales)
})

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

ETA_THRESH  <- 0.01   # 1% variance explained
N_TOP_LABEL <- 20

INPUT_CSV <- file.path(RESULTS_DIR, "variance_explained_probe_barcode_NTC.csv")

if (!file.exists(INPUT_CSV)) {
  stop("Eta-squared CSV not found: ", INPUT_CSV,
       "\nRun 04_compute_variance_explained.R first.")
}

# ==============================================================================
# Load and label
# ==============================================================================

df_eta <- read.csv(INPUT_CSV) |> arrange(rank_bc)

n_above <- sum(df_eta$eta_sq_bc >= ETA_THRESH)
message(sprintf("Genes with eta\u00b2 >= 1%%: %d / %d", n_above, nrow(df_eta)))

above_thresh  <- df_eta[df_eta$eta_sq_bc >= ETA_THRESH, ]
label_genes   <- head(above_thresh$gene, N_TOP_LABEL)
df_eta$label  <- ifelse(df_eta$gene %in% label_genes, df_eta$gene, NA_character_)

# ==============================================================================
# Plot
# ==============================================================================

p <- ggplot(df_eta, aes(x = rank_bc, y = eta_sq_bc,
                         color = eta_sq_bc >= ETA_THRESH)) +
  geom_point(alpha = 0.6, size = 1.2, shape = 16) +
  scale_color_manual(values = c("TRUE" = "tomato3", "FALSE" = "grey60")) +
  geom_hline(yintercept = ETA_THRESH,
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(
    aes(label = label),
    size = 7, force = 0.5, max.overlaps = 12,
    box.padding = 0.2, point.padding = 0.15, segment.size = 0.2,
    na.rm = TRUE, color = "black", segment.color = "black"
  ) +
  annotate("text",
    x = nrow(df_eta), y = ETA_THRESH, vjust = -0.5, hjust = 1,
    label = sprintf("\u03b7\u00b2 = 1%% (n = %d genes)", n_above),
    size = 7, color = "black") +
  scale_x_continuous(limits = c(1, nrow(df_eta))) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(df_eta$eta_sq_bc) * 1.05)) +
  labs(x = "Gene rank (by probe barcode \u03b7\u00b2)",
       y = "Variance explained by probe barcode") +
  theme_few() +
  theme(
    axis.title        = element_text(size = 20, color = "black"),
    axis.text         = element_text(size = 16, color = "black"),
    legend.position   = "none",
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "fig1f_ranked_eta_sq_probe_barcode.pdf"),
       p, width = 7, height = 5.8, device = cairo_pdf)
message("Saved: fig1f_ranked_eta_sq_probe_barcode.pdf")

message("Done.")
print(sessionInfo())
