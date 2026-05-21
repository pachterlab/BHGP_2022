library(tidyverse)

dat <- tibble(
  method      = factor(c("logp1\n(no Harmony)", "logp1\n+ Harmony", "CLR\n(no Harmony)", "CLR\n+ Harmony"),
                       levels = c("logp1\n(no Harmony)", "logp1\n+ Harmony",
                                  "CLR\n(no Harmony)",   "CLR\n+ Harmony")),
  knn_overlap = c(0.3751, 0.4299, 0.3874, 0.4790),
  pc1_depth_r = c(0.921,  0.824,  0.905,  0.489),
  norm        = c("logp1","logp1","CLR","CLR"),
  integration = c("none","Harmony","none","Harmony")
)

colors <- c(
  "logp1\n(no Harmony)" = "#6baed6",
  "logp1\n+ Harmony"    = "#2171b5",
  "CLR\n(no Harmony)"   = "#74c476",
  "CLR\n+ Harmony"      = "#238b45"
)

p_overlap <- ggplot(dat, aes(x = method, y = knn_overlap, fill = method)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.6), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL,
       y = "KNN overlap with\nfull-depth ground truth (k=50)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

p_confound <- ggplot(dat, aes(x = method, y = pc1_depth_r, fill = method)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1.02), expand = expansion(mult = c(0, 0.05)),
                     labels = scales::number_format(accuracy = 0.1)) +
  labs(x = NULL,
       y = "|r| PC1 ~ log(depth)\n(lower = less depth confounding)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

p_combined <- cowplot::plot_grid(
  p_overlap, p_confound,
  nrow = 1, align = "h", axis = "tb",
  labels = c("A", "B"), label_size = 12
)

title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "2,000 reads/cell downsampling: CLR + Harmony outperforms logp1 + Harmony\nmcSCRB: 50% full-depth + 50% downsampled cells, integrated with Harmony",
    fontface = "plain", size = 9.5, x = 0.5, hjust = 0.5
  )

final <- cowplot::plot_grid(title, p_combined, ncol = 1, rel_heights = c(0.12, 1))

ggsave("../output/clr_harmony_2k.pdf", final, width = 6, height = 3.8)
cat("Saved: ../output/clr_harmony_2k.pdf\n")
