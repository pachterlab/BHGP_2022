library(tidyverse)

results <- tribble(
  ~target_depth, ~method,               ~norm,   ~integration, ~knn_overlap, ~pc1_depth_r,
  500,  "logp1 (no Seurat)",   "logp1", "none",    0.3451,  0.980,
  500,  "logp1 + Seurat",      "logp1", "Seurat",  0.3916,  0.099,
  500,  "CLR (no Seurat)",     "CLR",   "none",    0.3433,  0.984,
  500,  "CLR + Seurat",        "CLR",   "Seurat",  0.3840,  0.107,
  1000, "logp1 (no Seurat)",   "logp1", "none",    0.3543,  0.969,
  1000, "logp1 + Seurat",      "logp1", "Seurat",  0.4238,  0.079,
  1000, "CLR (no Seurat)",     "CLR",   "none",    0.3544,  0.974,
  1000, "CLR + Seurat",        "CLR",   "Seurat",  0.4206,  0.065,
  2000, "logp1 (no Seurat)",   "logp1", "none",    0.3751,  0.921,
  2000, "logp1 + Seurat",      "logp1", "Seurat",  0.4729,  0.061,
  2000, "CLR (no Seurat)",     "CLR",   "none",    0.3874,  0.905,
  2000, "CLR + Seurat",        "CLR",   "Seurat",  0.4837,  0.043,
  5000, "logp1 (no Seurat)",   "logp1", "none",    0.4708,  0.537,
  5000, "logp1 + Seurat",      "logp1", "Seurat",  0.4994,  0.088,
  5000, "CLR (no Seurat)",     "CLR",   "none",    0.5160,  0.184,
  5000, "CLR + Seurat",        "CLR",   "Seurat",  0.5039,  0.033,
  10000,"logp1 (no Seurat)",   "logp1", "none",    0.6173,  0.166,
  10000,"logp1 + Seurat",      "logp1", "Seurat",  0.5214,  0.099,
  10000,"CLR (no Seurat)",     "CLR",   "none",    0.6504,  0.010,
  10000,"CLR + Seurat",        "CLR",   "Seurat",  0.5343,  0.041,
) %>%
  mutate(
    method = factor(method, levels = c("logp1 (no Seurat)", "logp1 + Seurat",
                                        "CLR (no Seurat)",   "CLR + Seurat")),
    pc1_depth_r = abs(pc1_depth_r)
  )

colors <- c(
  "logp1 (no Seurat)" = "#6baed6",
  "logp1 + Seurat"    = "#2171b5",
  "CLR (no Seurat)"   = "#74c476",
  "CLR + Seurat"      = "#238b45"
)
shapes <- c(
  "logp1 (no Seurat)" = 1,
  "logp1 + Seurat"    = 16,
  "CLR (no Seurat)"   = 0,
  "CLR + Seurat"      = 15
)
linetypes <- c(
  "logp1 (no Seurat)" = "dashed",
  "logp1 + Seurat"    = "solid",
  "CLR (no Seurat)"   = "dashed",
  "CLR + Seurat"      = "solid"
)

p_overlap <- ggplot(results, aes(x = target_depth, y = knn_overlap,
                                  color = method, shape = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = c(500, 1000, 2000, 5000, 10000),
                labels = c("500", "1k", "2k", "5k", "10k")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.3, 0.72)) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = linetypes) +
  labs(x = "Downsampled depth (reads/cell)",
       y = "KNN overlap with\nfull-depth ground truth (k=50)",
       color = NULL, shape = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

p_confound <- ggplot(results, aes(x = target_depth, y = pc1_depth_r,
                                   color = method, shape = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = c(500, 1000, 2000, 5000, 10000),
                labels = c("500", "1k", "2k", "5k", "10k")) +
  scale_y_continuous(limits = c(0, 1.02),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = linetypes) +
  labs(x = "Downsampled depth (reads/cell)",
       y = "|r| PC1 ~ log(depth)\n(lower = less depth confounding)",
       color = NULL, shape = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 2))

p_combined <- cowplot::plot_grid(
  p_overlap, p_confound,
  nrow = 2, align = "v", axis = "lr",
  labels = c("A", "B"), label_size = 12
)

title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "CLR + Seurat integration vs logp1 + Seurat\nmcSCRB: 50% full-depth + 50% downsampled cells, integrated with Seurat CCA",
    fontface = "plain", size = 10, x = 0.5, hjust = 0.5
  )

final <- cowplot::plot_grid(title, p_combined, ncol = 1, rel_heights = c(0.08, 1))

ggsave("../output/clr_seurat_sweep.pdf", final, width = 4.5, height = 6)
cat("Saved: ../output/clr_seurat_sweep.pdf\n")
