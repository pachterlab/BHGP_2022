#!/usr/bin/env Rscript
# plot_fig2d_clr.R
#
# Reproduces Figure 2d (siRNA KD downsampling panel) with PFlog1pPF (CLR) added.
# Run from notebooks/ directory: Rscript plot_fig2d_clr.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

source("utils.R")
source("annotation_helper.R")

# Rename CLR to PFlog1pPF (CLR)
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

# ── Build pca_dat (same pipeline as published notebook) ──────────────────────
set.seed(42)  # for fct_reorder(runif)
pca_dat <- res %>%
  filter(alpha %in% c("TRUE", "FALSE")) %>%
  filter(knn == 50) %>%
  mutate(knn_recovery = overlap) %>%
  group_by(benchmark, dataset, pca_dim, transformation) %>%
  summarize(knn_recovery = mean(knn_recovery), .groups = "drop") %>%
  left_join(trans_families, by = "transformation") %>%
  mutate(transformation = fct_reorder(transformation, runif(n()))) %>%
  mutate(is_dim_indep = transformation == "sanity_dists",
         is_clr       = transformation == "clr")

# siRNA KD data
sirna_dat <- pca_dat %>% filter(dataset == "smartSeq3_siRNA_knockdown")

# Endpoint labels for highlighted lines
label_pts <- sirna_dat %>%
  filter(transformation %in% c("clr", "scgpt")) %>%
  group_by(transformation) %>%
  slice_max(pca_dim, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = c(clr = "PFlog1pPF (CLR)", scgpt = "scGPT")[as.character(transformation)])

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(sirna_dat, aes(x = pca_dim, y = knn_recovery)) +
  # All non-CLR lines (thin)
  geom_line(data = filter(sirna_dat, !is_clr),
            aes(group = transformation, color = family, linetype = is_dim_indep),
            linewidth = 0.4, show.legend = FALSE) +
  # CLR line (thicker, same color family = teal)
  geom_line(data = filter(sirna_dat, is_clr),
            aes(group = transformation, color = family),
            linewidth = 1.2, show.legend = FALSE) +
  # Endpoint labels
  geom_text(data = label_pts,
            aes(label = label, color = family),
            hjust = -0.08, vjust = 0.5,
            size = font_size_small / ggplot2::.pt,
            show.legend = FALSE) +
  scale_x_log10(breaks = c(5, 10, 50, 100)) +
  scale_y_continuous(limits = c(0, 50),
                     breaks = c(0, 10, 20, 30, 40, 50)) +
  scale_color_manual(values = trans_families_colors) +
  scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "solid")) +
  coord_cartesian(clip = "off") +
  labs(
    title = "siRNA KD (Downsampling)",
    y     = "k-NN\nOverlap",
    x     = "No. PCA-dimensions"
  ) +
  theme(
    plot.title       = element_text(face = "plain", size = font_size_small,
                                    hjust = 0.5, margin = margin(), lineheight = 0.8),
    axis.title       = element_text(size = font_size_small),
    axis.text        = element_text(size = font_size_small),
    plot.margin      = margin(5, 60, 5, 5)
  )

ggsave("fig2d_siRNA_clr.png", p, width = 4, height = 3, dpi = 300, bg = "white")
message("Saved: fig2d_siRNA_clr.png")
