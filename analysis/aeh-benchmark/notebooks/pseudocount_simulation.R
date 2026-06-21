#!/usr/bin/env Rscript
# Supplementary pseudocount simulation.
#
# The old version of this figure plotted an abstract "correction factor" cf,
# where the shifted-log pseudocount was 1 / (cf * alpha).  This version plots the
# same idea in the notation used in the paper:
#
#   log(1 + K x_gc / s_c)
#
# where s_c is the cell depth.  On the count scale the corresponding pseudocount
# is approximately sbar / K, so the delta-method / Anscombe choice
# y0 = 1 / (4 alpha) corresponds to
#
#   K = 4 alpha sbar.
#
# The simulation below generates Huber-benchmark-like settings with fixed
# underlying proportions and variable cell depths.  For a range of pseudocount
# values, equivalently for a range of K values with y0 = sbar / K, it
# evaluates the modeled delta-method technical variance
#
#   q_gc(K) = {d/dmu log(1 + K mu_gc / s_c)}^2 (mu_gc + alpha mu_gc^2)
#
# and summarizes how evenly q_g(K) = mean_c q_gc(K) is stabilized across genes.
# Lower CV is better.  The plotted minimum is the minimum of this finite
# gene/depth distribution; the Anscombe line is the analytic shifted-log choice.
#
# Usage: Rscript pseudocount_simulation.R
# Writes: ../output/pseudocount_simulation.{pdf,png,csv}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

set.seed(1)

n_genes <- 2000L
n_cells <- 500L
n_rep <- 20L
angelidis_alpha <- 0.4047020313316354

scenarios <- tibble::tribble(
  ~scenario_id, ~alpha, ~mean_depth, ~alpha_label,
  "ss3_fib", 0.0903663957233413, 228770.623306233, "0.090 (Smart-seq3 fibroblasts)",
  "gsm4235299", 0.23372025889020434, 55641.95855614973, "0.234 (GSM4235299)",
  "sirna", 0.242038898430308, 125183.563983249, "0.242 (siRNA KD)",
  "angelidis", angelidis_alpha, 818.46, "0.405 (Angelidis et al.)",
  "mcscrb", 0.796583324572214, 63974.1124497992, "0.797 (mcSCRB)",
  "ss3_hek", 1.06442469733767, 50291.0427350427, "1.06 (Smart-seq3 HEK)",
  "gsm5429730", 4.030432793492418, 7401.56266768877, "4.03 (GSM5429730)",
  "gsm5429729", 6.64734493432392, 7254.195482034863, "6.65 (GSM5429729)"
)

k_reference <- c(
  "CP10k" = 1e4,
  "CPM" = 1e6
)

make_k_values <- function(alpha, mean_depth) {
  sort(unique(c(
    10^seq(2, 7, length.out = 81),
    unname(k_reference),
    4 * alpha * mean_depth
  )))
}

cv_delta_q <- function(gene_prop, depth, alpha, K) {
  # mu_gc = p_g s_c, so after depth normalization mu_gc / s_c = p_g.
  # q_gc varies across cells because the count-scale NB variance still depends
  # on s_c.
  denom <- 1 + K * gene_prop
  slope2 <- outer((K / depth)^2, 1 / denom^2)
  mu <- outer(depth, gene_prop)
  q <- slope2 * (mu + alpha * mu^2)
  q_gene <- colMeans(q)
  sd(q_gene) / mean(q_gene)
}

simulate_once <- function(alpha, mean_depth) {
  gene_weight <- 10^runif(n_genes, -1, 1)
  gene_prop <- gene_weight / sum(gene_weight)

  depth <- round(mean_depth * exp(rnorm(n_cells, sd = 0.25)))
  depth <- pmax(depth, 100L)

  k_values <- make_k_values(alpha, mean_depth)

  data.frame(
    K = k_values,
    cv = vapply(k_values, function(K) cv_delta_q(gene_prop, depth, alpha, K), numeric(1))
  )
}

results <- bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
  scenario <- scenarios[i, ]
  bind_rows(lapply(seq_len(n_rep), function(rep_id) {
    simulate_once(scenario$alpha, scenario$mean_depth) |>
      mutate(
        scenario_id = scenario$scenario_id,
        alpha = scenario$alpha,
        mean_depth = scenario$mean_depth,
        alpha_label = scenario$alpha_label,
        replicate = rep_id
      )
  }))
}))

summary_df <- results |>
  group_by(scenario_id, alpha, mean_depth, alpha_label, K) |>
  summarize(
    median_cv = median(cv),
    q25 = quantile(cv, 0.25),
    q75 = quantile(cv, 0.75),
    .groups = "drop"
  ) |>
  mutate(
    pseudocount = mean_depth / K,
    alpha_label = factor(alpha_label, levels = scenarios$alpha_label)
  ) |>
  arrange(alpha_label, pseudocount)

anscombe_df <- scenarios |>
  transmute(
    scenario_id,
    alpha_label = factor(alpha_label, levels = scenarios$alpha_label),
    pseudocount = 1 / (4 * alpha),
    scenario = "Anscombe"
  )

min_df <- summary_df |>
  slice_min(median_cv, by = scenario_id, n = 1, with_ties = FALSE) |>
  transmute(
    scenario_id,
    alpha_label,
    pseudocount,
    scenario = "CV minimum"
  )

ref_df <- merge(
  scenarios,
  data.frame(scenario = names(k_reference), K = unname(k_reference)),
  by = NULL
) |>
  transmute(
    scenario_id,
    alpha_label = factor(alpha_label, levels = scenarios$alpha_label),
    pseudocount = mean_depth / K,
    scenario
  )

out <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "..", "output")
if (!dir.exists(out)) out <- "../output"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

write.csv(results, file.path(out, "pseudocount_simulation.csv"), row.names = FALSE)

line_cols <- c(
  "Anscombe" = "#aa3377",
  "CV minimum" = "#d55e00",
  "CP10k" = "#333333",
  "CPM" = "#666666"
)

p <- ggplot(summary_df, aes(x = pseudocount, y = median_cv)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#b9c9d6", alpha = 0.55) +
  geom_line(color = "#2f5d7c", linewidth = 0.9) +
  geom_vline(
    data = ref_df,
    aes(xintercept = pseudocount, color = scenario),
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  geom_vline(
    data = anscombe_df,
    aes(xintercept = pseudocount, color = scenario),
    linewidth = 0.65
  ) +
  geom_vline(
    data = min_df,
    aes(xintercept = pseudocount, color = scenario),
    linewidth = 0.75
  ) +
  facet_wrap(vars(alpha_label), nrow = 2) +
  scale_x_log10(
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
    labels = c(expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), "1", "10", "100", "1000")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_color_manual(
    values = line_cols,
    breaks = c("Anscombe", "CV minimum", "CP10k", "CPM"),
    name = NULL
  ) +
  labs(
    x = "count-scale pseudocount",
    y = expression("CV of modeled technical variances " * q)
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f2f2f2", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(size = 9)
  ) +
  guides(color = guide_legend(override.aes = list(
    linewidth = c(0.65, 0.75, 0.45, 0.45),
    linetype = c("solid", "solid", "dashed", "dashed")
  )))

ggsave(file.path(out, "pseudocount_simulation.pdf"), p, width = 10.2, height = 6.1)
ggsave(file.path(out, "pseudocount_simulation.png"), p, width = 10.2, height = 6.1, dpi = 300)

cat("wrote", file.path(out, "pseudocount_simulation.{pdf,png,csv}"), "\n")
