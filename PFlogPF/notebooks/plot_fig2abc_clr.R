#!/usr/bin/env Rscript
# plot_fig2abc_clr.R
#
# Reproduces Figure 2 panels a, b, c from Ahlmann-Eltze & Huber (2023)
# with CLR (centered log-ratio) added as an additional transformation.
#
# Run from notebooks/ directory:
#   Rscript plot_fig2abc_clr.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggbeeswarm)
})

# ── Stubs for gggroupedscale (not available) ──────────────────────────────────
scale_y_grouped_discrete <- function(..., grouping = NULL, gap_size = 1, limits = NULL,
                                     add_group_label = FALSE, guide = NULL) {
  args <- list(...)
  if (!is.null(limits)) args[["limits"]] <- limits
  if (!is.null(guide)) args[["guide"]] <- guide_axis()
  do.call(scale_y_discrete, args)
}
theme_grouped_axis  <- function(...) theme()
guide_grouped_axis  <- function(...) guide_axis(...)

source("utils.R")
source("annotation_helper.R")

# Rename CLR to PFlog1pPF (CLR) in plain labels
trans_labels_plain["clr"] <- "PFlog1pPF (CLR)"

# ── Load data ─────────────────────────────────────────────────────────────────
res <- bind_rows(
  read_tsv("../benchmark/output/benchmark_results/simulation_results.tsv",
           show_col_types = FALSE) %>%
    transmute(benchmark = "simulation", overlap = mean_knn_overlap, knn, pca_dim,
              alpha = as.character(alpha), transformation, dataset = simulator, replicate = seed),
  read_tsv("../benchmark/output/benchmark_results/consistency_results.tsv",
           show_col_types = FALSE) %>%
    transmute(benchmark = "consistency", overlap = mean_overlap, knn, pca_dim,
              alpha = as.character(alpha), transformation, dataset, replicate = seed),
  read_tsv("../benchmark/output/benchmark_results/downsampling_results.tsv",
           show_col_types = FALSE) %>%
    transmute(benchmark = "downsampling", overlap = overlap, knn, pca_dim,
              alpha = as.character(alpha), transformation, dataset, replicate = seed)
) %>%
  mutate(transformation = factor(transformation, levels = trans_families$transformation))

# Same parameter choices as the published figure
parameter_choices <- bind_rows(
  tibble(benchmark = "downsampling", knn = 50, pca_dim = c(10, 10, 10, 10, 50),
         dataset = c("mcSCRB", "smartSeq3_fibroblasts", "smartSeq3_fibroblasts_alt",
                     "smartSeq3_hek", "smartSeq3_siRNA_knockdown")),
  tibble(benchmark = "simulation", knn = 50, pca_dim = c(5, 10, 10, 200, 50),
         dataset = c("dyngen", "linear_walk", "muscat", "random_walk", "scDesign2")),
  tibble(benchmark = "consistency", knn = 50, pca_dim = 50,
         dataset = unique(filter(res, benchmark == "consistency")$dataset))
)

# ── Build res_main ────────────────────────────────────────────────────────────
res_main <- res %>%
  filter(alpha %in% c("TRUE", "FALSE")) %>%
  inner_join(parameter_choices, by = c("benchmark", "knn", "pca_dim", "dataset")) %>%
  mutate(knn_recovery = overlap / knn) %>%
  group_by(dataset, replicate, knn) %>%
  mutate(knn_recovery = knn_recovery / mean(knn_recovery)) %>%
  ungroup() %>%
  left_join(trans_families, by = "transformation")

# ── Per-panel axis limits ─────────────────────────────────────────────────────
x_ranges <- res_main %>%
  group_by(benchmark) %>%
  summarize(xmax = max(knn_recovery), xmin = min(knn_recovery), .groups = "drop")

get_xlim <- function(bm) {
  xmax <- x_ranges$xmax[x_ranges$benchmark == bm]
  xmin <- x_ranges$xmin[x_ranges$benchmark == bm]
  c(min(0.0, floor(xmin * 10) / 10), ceiling(xmax * 1.05 * 2) / 2)
}

make_breaks <- function(xlim) {
  span <- diff(xlim)
  if (span <= 2)      seq(0, xlim[2], by = 0.5)
  else if (span <= 5) seq(0, xlim[2], by = 1)
  else                seq(0, xlim[2], by = 2)
}

# ── Family bracket helper ─────────────────────────────────────────────────────
family_bracket_data <- function(data) {
  levs <- rev(levels(droplevels(factor(data$transformation,
                                       levels = trans_families$transformation))))
  tibble(transformation = levs, ypos = seq_len(length(levs))) %>%
    left_join(trans_families, by = "transformation") %>%
    group_by(family) %>%
    summarize(ymin = min(ypos) - 0.45, ymax = max(ypos) + 0.45,
              ymid = (min(ypos) + max(ypos)) / 2, .groups = "drop") %>%
    mutate(
      color = trans_families_colors[family],
      label = as.character(trans_families_labels[family])
    )
}

# ── Plot function ─────────────────────────────────────────────────────────────
make_panel <- function(data, bm, add_group_label = FALSE, subtitle = "") {
  xlim   <- get_xlim(bm)
  breaks <- make_breaks(xlim)
  xspan  <- diff(xlim)
  brack_x <- xlim[1] - xspan * 0.70
  label_x <- xlim[1] - xspan * 0.84

  p <- ggplot(data, aes(x = knn_recovery, y = transformation, color = family)) +
    geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, color = "grey40") +
    ggbeeswarm::geom_quasirandom(color = "grey75", size = 0.9, alpha = 0.6,
                                 groupOnX = FALSE) +
    stat_summary(geom = "point", position = position_dodge2(width = 0.35),
                 fun = mean, size = 2.2) +
    scale_y_grouped_discrete(
      grouping        = ~ trans_families_labels[deframe(trans_families)[.x]],
      gap_size        = 1.3,
      limits          = rev,
      labels          = trans_labels_plain,
      add_group_label = add_group_label,
      guide           = if (add_group_label) guide_grouped_axis() else guide_axis()
    ) +
    scale_x_continuous(breaks = breaks) +
    coord_cartesian(xlim = xlim, clip = "off") +
    scale_color_manual(values = trans_families_colors, labels = trans_families_labels,
                       guide = "none") +
    theme_grouped_axis() +
    theme(
      axis.title.y        = element_blank(),
      plot.title.position = "plot",
      plot.subtitle       = element_text(size = font_size_small),
      axis.text.y         = element_text(size = font_size_small),
      axis.text.x         = element_text(size = font_size_small),
      axis.title.x        = element_text(size = font_size_small),
      plot.margin         = if (add_group_label) margin(5, 5, 5, 200) else margin(5, 5, 5, 5),
      panel.grid.major.x  = element_line(color = "grey92", linewidth = 0.3)
    ) +
    labs(x = "Relative k-NN Overlap", subtitle = subtitle)

  if (add_group_label) {
    fd <- family_bracket_data(data)
    for (i in seq_len(nrow(fd))) {
      p <- p +
        annotate("segment", x = brack_x, xend = brack_x,
                 y = fd$ymin[i], yend = fd$ymax[i],
                 color = fd$color[i], linewidth = 1.3) +
        annotate("text", x = label_x, y = fd$ymid[i],
                 label = fd$label[i], angle = 90,
                 size = font_size_tiny / ggplot2::.pt,
                 color = fd$color[i], hjust = 0.5, vjust = 0.5)
    }
  }
  p
}

# ── Build panels ──────────────────────────────────────────────────────────────
panel_a <- res_main %>% filter(benchmark == "consistency") %>%
  make_panel("consistency", add_group_label = TRUE,
             subtitle = "10X gene subset 1 vs. gene subset 2")

panel_b <- res_main %>% filter(benchmark == "simulation") %>%
  make_panel("simulation", subtitle = "Ground truth vs. simulated counts")

panel_c <- res_main %>% filter(benchmark == "downsampling") %>%
  make_panel("downsampling", subtitle = "Original vs. downsampled deep-seq data")

# ── Combine and save ──────────────────────────────────────────────────────────
fig <- plot_grid(panel_a, panel_b, panel_c,
                 nrow = 1, labels = NULL,
                 rel_widths = c(1.6, 1, 1))

ggsave("fig2abc_with_clr.png", fig, width = 11, height = 7.5, dpi = 300, bg = "white")
message("Saved: fig2abc_with_clr.png")
