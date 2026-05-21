library(tidyverse)

results <- tribble(
  ~target_depth, ~method,                  ~norm,   ~integration, ~knn_overlap, ~pc1_depth_r,
  500,  "logp1 (no Harmony)",  "logp1", "none",    0.3333, 0.549,
  500,  "logp1 + Harmony",     "logp1", "Harmony", 0.4377, 0.275,
  500,  "CLR (no Harmony)",    "CLR",   "none",    0.3322, 0.407,
  500,  "CLR + Harmony",       "CLR",   "Harmony", 0.4337, 0.253,
  1000, "logp1 (no Harmony)",  "logp1", "none",    0.3859, 0.423,
  1000, "logp1 + Harmony",     "logp1", "Harmony", 0.5327, 0.310,
  1000, "CLR (no Harmony)",    "CLR",   "none",    0.3900, 0.369,
  1000, "CLR + Harmony",       "CLR",   "Harmony", 0.5275, 0.298,
  2000, "logp1 (no Harmony)",  "logp1", "none",    0.4618, 0.447,
  2000, "logp1 + Harmony",     "logp1", "Harmony", 0.6404, 0.395,
  2000, "CLR (no Harmony)",    "CLR",   "none",    0.4741, 0.418,
  2000, "CLR + Harmony",       "CLR",   "Harmony", 0.6360, 0.384,
  5000, "logp1 (no Harmony)",  "logp1", "none",    0.6760, 0.565,
  5000, "logp1 + Harmony",     "logp1", "Harmony", 0.8136, 0.555,
  5000, "CLR (no Harmony)",    "CLR",   "none",    0.6993, 0.548,
  5000, "CLR + Harmony",       "CLR",   "Harmony", 0.8008, 0.542,
  7000, "logp1 (no Harmony)",  "logp1", "none",    0.8246, 0.605,
  7000, "logp1 + Harmony",     "logp1", "Harmony", 0.8822, 0.603,
  7000, "CLR (no Harmony)",    "CLR",   "none",    0.8286, 0.590,
  7000, "CLR + Harmony",       "CLR",   "Harmony", 0.8621, 0.589,
) %>%
  mutate(
    method = factor(method, levels = c("logp1 (no Harmony)", "logp1 + Harmony",
                                        "CLR (no Harmony)",   "CLR + Harmony")),
    pc1_depth_r = abs(pc1_depth_r)
  )

colors <- c(
  "logp1 (no Harmony)" = "#6baed6",
  "logp1 + Harmony"    = "#2171b5",
  "CLR (no Harmony)"   = "#74c476",
  "CLR + Harmony"      = "#238b45"
)
shapes <- c(
  "logp1 (no Harmony)" = 1,
  "logp1 + Harmony"    = 16,
  "CLR (no Harmony)"   = 0,
  "CLR + Harmony"      = 15
)
linetypes <- c(
  "logp1 (no Harmony)" = "dashed",
  "logp1 + Harmony"    = "solid",
  "CLR (no Harmony)"   = "dashed",
  "CLR + Harmony"      = "solid"
)

p_overlap <- ggplot(results, aes(x = target_depth, y = knn_overlap,
                                  color = method, shape = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = c(500, 1000, 2000, 5000, 7000),
                labels = c("500", "1k", "2k", "5k", "7k")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.3, 0.95)) +
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
  scale_x_log10(breaks = c(500, 1000, 2000, 5000, 7000),
                labels = c("500", "1k", "2k", "5k", "7k")) +
  scale_y_continuous(limits = c(0, 0.65),
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
    "CLR + Harmony — Mouse Aorta 10x (GSE178765, ~9,400 reads/cell full depth)\n50% full-depth + 50% downsampled cells, integrated with Harmony",
    fontface = "plain", size = 10, x = 0.5, hjust = 0.5
  )

final <- cowplot::plot_grid(title, p_combined, ncol = 1, rel_heights = c(0.08, 1))

ggsave("../output/clr_harmony_10x_sweep.pdf", final, width = 4.5, height = 6)
cat("Saved: ../output/clr_harmony_10x_sweep.pdf\n")

# ── 2k barplot ────────────────────────────────────────────────────────────────
dat2k <- results %>% filter(target_depth == 2000) %>%
  mutate(method_lab = factor(
    c("logp1\n(no Harmony)", "logp1\n+ Harmony", "CLR\n(no Harmony)", "CLR\n+ Harmony"),
    levels = c("logp1\n(no Harmony)", "logp1\n+ Harmony", "CLR\n(no Harmony)", "CLR\n+ Harmony")
  ))

colors2 <- setNames(unname(colors), levels(dat2k$method_lab))

p_ov2k <- ggplot(dat2k, aes(x = method_lab, y = knn_overlap, fill = method_lab)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = colors2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.75), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "KNN overlap with\nfull-depth ground truth (k=50)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

p_r2k <- ggplot(dat2k, aes(x = method_lab, y = pc1_depth_r, fill = method_lab)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = colors2) +
  scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.05)),
                     labels = scales::number_format(accuracy = 0.1)) +
  labs(x = NULL, y = "|r| PC1 ~ log(depth)\n(lower = less depth confounding)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

p2k <- cowplot::plot_grid(p_ov2k, p_r2k, nrow = 1, align = "h", axis = "tb",
                           labels = c("A", "B"), label_size = 12)
title2k <- cowplot::ggdraw() +
  cowplot::draw_label(
    "2,000 reads/cell — Mouse Aorta 10x (GSE178765)\n50% full-depth + 50% downsampled cells, integrated with Harmony",
    fontface = "plain", size = 9.5, x = 0.5, hjust = 0.5
  )
final2k <- cowplot::plot_grid(title2k, p2k, ncol = 1, rel_heights = c(0.12, 1))
ggsave("../output/clr_harmony_10x_2k.pdf", final2k, width = 6, height = 3.8)
cat("Saved: ../output/clr_harmony_10x_2k.pdf\n")
