#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggh4x)
  library(cowplot)
  library(rlang)
})

script_path <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]])))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  stop("Could not determine script path", call. = FALSE)
}

script_file <- script_path()
repo_dir <- normalizePath(file.path(dirname(script_file), ".."))
pf_dir <- file.path(repo_dir, "analysis", "aeh-benchmark")
notebook_dir <- file.path(pf_dir, "notebooks")
benchmark_dir <- file.path(pf_dir, "benchmark", "output", "benchmark_results")
figure_dir <- normalizePath(file.path(repo_dir, "..", "figures"), mustWork = FALSE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Compatibility shim used by the benchmark plotting notebook.
# gggroupedscale is not available for current ggplot2, so use ggh4x nested axes.
scale_y_grouped_discrete <- function(..., grouping = NULL, gap_size = 1, limits = NULL,
                                     add_group_label = FALSE, guide = NULL,
                                     labels = NULL) {
  args <- list(...)
  if (!is.null(limits)) args[["limits"]] <- limits
  if (add_group_label && !is.null(grouping) && !is.null(labels)) {
    methods <- names(labels)
    grouping_fn <- if (inherits(grouping, "formula")) {
      rlang::as_function(grouping)
    } else if (is.function(grouping)) {
      grouping
    } else {
      function(.x) grouping[.x]
    }
    fams <- vapply(methods, function(m) as.character(grouping_fn(m)), character(1))
    args[["labels"]] <- setNames(paste0(fams, " | ", labels), methods)
    args[["guide"]] <- ggh4x::guide_axis_nested(delim = "|", n.dodge = 1)
  } else {
    args[["labels"]] <- labels
    args[["guide"]] <- guide_axis()
  }
  do.call(scale_y_discrete, args)
}

theme_grouped_axis <- function(...) theme()
guide_grouped_axis <- function(...) guide_axis(...)
options(ggrastr.default.dpi = 150)

source(file.path(notebook_dir, "utils.R"), chdir = TRUE)
source(file.path(notebook_dir, "annotation_helper.R"), chdir = TRUE)

# Keep this adapted Supplementary Note figure visually aligned with the
# Ahlmann-Eltze--Huber benchmark style.
trans_labels[c("logp1", "logp1_hvg", "logp1_hvg_zscore", "logp1_zscore", "logp1_size_normed")] <- c(
  "logp1" = r"($\log(y/s+1)$)",
  "logp1_hvg" = r"($\log(y/s+1)\rightarrow$HVG)",
  "logp1_hvg_zscore" = r"($\log(y/s+1)\rightarrow$HVG$\rightarrow$Z)",
  "logp1_zscore" = r"($\log(y/s+1)\rightarrow$Z)",
  "logp1_size_normed" = r"($\log(y/s+1)/u$)"
)

read_benchmark_results <- function() {
  bind_rows(
    readr::read_tsv(file.path(benchmark_dir, "simulation_results.tsv"), show_col_types = FALSE) %>%
      transmute(
        benchmark = "simulation",
        overlap = mean_knn_overlap,
        knn,
        pca_dim,
        alpha = as.character(alpha),
        transformation,
        dataset = simulator,
        replicate = seed,
        cputime_sec,
        elapsed_sec
      ),
    readr::read_tsv(file.path(benchmark_dir, "consistency_results.tsv"), show_col_types = FALSE) %>%
      transmute(
        benchmark = "consistency",
        overlap = mean_overlap,
        knn,
        pca_dim,
        alpha = as.character(alpha),
        transformation,
        dataset,
        replicate = seed,
        cputime_sec,
        elapsed_sec
      ),
    readr::read_tsv(file.path(benchmark_dir, "downsampling_results.tsv"), show_col_types = FALSE) %>%
      transmute(
        benchmark = "downsampling",
        overlap,
        knn,
        pca_dim,
        alpha = as.character(alpha),
        transformation,
        dataset,
        replicate = seed,
        cputime_sec = full_cputime_sec,
        elapsed_sec = full_elapsed_sec
      )
  ) %>%
    mutate(transformation = factor(transformation, levels = trans_families$transformation))
}

read_clr_overrides <- function() {
  override_files <- file.path(
    figure_dir,
    c(
      "fig2abc_clr_alpha_k_simulation.tsv",
      "fig2abc_clr_alpha_k_consistency.tsv",
      "fig2abc_clr_alpha_k_downsampling.tsv"
    )
  )
  missing_files <- override_files[!file.exists(override_files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing scclr CLR override table(s):\n",
      paste(missing_files, collapse = "\n"),
      "\nRun scripts/recompute_fig2abc_scclr.py first.",
      call. = FALSE
    )
  }

  bind_rows(
    readr::read_tsv(file.path(figure_dir, "fig2abc_clr_alpha_k_simulation.tsv"), show_col_types = FALSE) %>%
      transmute(
        benchmark = "simulation",
        overlap = mean_knn_overlap,
        knn,
        pca_dim,
        alpha = as.character(alpha),
        transformation,
        dataset = simulator,
        replicate = seed,
        cputime_sec = NA_real_,
        elapsed_sec = NA_real_
      ),
    readr::read_tsv(file.path(figure_dir, "fig2abc_clr_alpha_k_consistency.tsv"), show_col_types = FALSE) %>%
      transmute(
        benchmark = "consistency",
        overlap = mean_overlap,
        knn,
        pca_dim,
        alpha = as.character(alpha),
        transformation,
        dataset,
        replicate = seed,
        cputime_sec = NA_real_,
        elapsed_sec = NA_real_
      ),
    readr::read_tsv(file.path(figure_dir, "fig2abc_clr_alpha_k_downsampling.tsv"), show_col_types = FALSE) %>%
      transmute(
        benchmark = "downsampling",
        overlap,
        knn,
        pca_dim,
        alpha = as.character(alpha),
        transformation,
        dataset,
        replicate = seed,
        cputime_sec = NA_real_,
        elapsed_sec = NA_real_
      )
  ) %>%
    mutate(transformation = factor(transformation, levels = trans_families$transformation))
}

replace_clr_rows <- function(res, overrides) {
  key_cols <- c("benchmark", "dataset", "replicate", "pca_dim", "knn")
  override_keys <- overrides %>%
    distinct(across(all_of(key_cols))) %>%
    mutate(.replace_clr = TRUE)

  res %>%
    left_join(override_keys, by = key_cols) %>%
    filter(!(transformation == "clr" & !is.na(.replace_clr))) %>%
    select(-.replace_clr) %>%
    bind_rows(overrides) %>%
    mutate(transformation = factor(transformation, levels = trans_families$transformation))
}

parameter_choices <- function(res) {
  bind_rows(
    tibble(
      benchmark = "downsampling",
      knn = 50,
      pca_dim = c(10, 10, 10, 10, 50),
      dataset = c(
        "mcSCRB",
        "smartSeq3_fibroblasts",
        "smartSeq3_fibroblasts_alt",
        "smartSeq3_hek",
        "smartSeq3_siRNA_knockdown"
      )
    ),
    tibble(
      benchmark = "simulation",
      knn = 50,
      pca_dim = c(5, 10, 50),
      dataset = c("dyngen", "muscat", "scDesign2")
    ),
    tibble(
      benchmark = "consistency",
      knn = 50,
      pca_dim = 50,
      dataset = unique(filter(res, benchmark == "consistency")$dataset)
    )
  )
}

make_main_plot_panel <- function(data, add_group_label = FALSE) {
  ggplot(data, aes(x = knn_recovery, y = transformation, color = family)) +
    geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2) +
    ggbeeswarm::geom_quasirandom(
      color = "grey",
      size = 0.4,
      alpha = 0.7,
      groupOnX = FALSE
    ) +
    stat_summary(geom = "point", fun = mean, size = 2.5, stroke = 0) +
    scale_y_grouped_discrete(
      grouping = ~ trans_families_labels[deframe(trans_families)[.x]],
      gap_size = 1.3,
      limits = rev,
      labels = trans_labels,
      add_group_label = add_group_label,
      guide = if (add_group_label) guide_grouped_axis() else guide_axis()
    ) +
    scale_x_continuous(breaks = c(0.5, 1, 1.5)) +
    coord_cartesian(xlim = c(0.2, 1.8)) +
    scale_color_manual(
      values = trans_families_colors,
      labels = trans_families_labels,
      guide = "none"
    ) +
    theme_grouped_axis(
      axis.grouping.line_padding = unit(5, "pt"),
      axis.grouping.line_height = unit(10, "pt"),
      axis.grouping.label.y = element_text(size = font_size_small, angle = 90)
    ) +
    theme(axis.title.y = element_blank(), plot.title.position = "plot")
}

res <- read_benchmark_results() %>%
  replace_clr_rows(read_clr_overrides())

choices <- parameter_choices(res)

res_main <- res %>%
  filter(alpha %in% c("TRUE", "FALSE")) %>%
  inner_join(choices, by = c("benchmark", "knn", "pca_dim", "dataset")) %>%
  mutate(knn_recovery = overlap / knn) %>%
  group_by(dataset, replicate, knn) %>%
  mutate(knn_recovery = knn_recovery / mean(knn_recovery)) %>%
  ungroup() %>%
  left_join(trans_families, by = "transformation")

clr_rows <- res_main %>% filter(transformation == "clr")
if (!identical(sort(unique(clr_rows$alpha)), "TRUE")) {
  stop("Expected all plotted CLR rows to use estimated alpha/K rows", call. = FALSE)
}

readr::write_tsv(
  res_main,
  file.path(figure_dir, "fig2abc_alpha_k_benchmark_style_source.tsv")
)

consistency_pl <- res_main %>%
  filter(benchmark == "consistency") %>%
  make_main_plot_panel(add_group_label = TRUE) +
  labs(
    x = "Relative $k$-NN overlap",
    subtitle = "10x gene subset 1 versus gene subset 2"
  )

simulation_pl <- res_main %>%
  filter(benchmark == "simulation") %>%
  make_main_plot_panel(add_group_label = FALSE) +
  coord_cartesian(xlim = c(0.2, 1.45)) +
  scale_x_continuous(breaks = c(0.45, 1, 1.25), labels = c("0.5", "1", "1.25")) +
  labs(
    x = "Relative $k$-NN overlap",
    subtitle = "Ground truth versus simulated counts"
  )

downsampling_pl <- res_main %>%
  filter(benchmark == "downsampling") %>%
  make_main_plot_panel(add_group_label = FALSE) +
  coord_cartesian(xlim = c(0.2, 3.8)) +
  scale_x_continuous(breaks = c(0.45, 1, 2, 3, 3.5), labels = c("0.5", "1", "2", "3", "3.5")) +
  labs(
    x = "Relative $k$-NN overlap",
    subtitle = "Original versus downsampled deep-seq data"
  )

out_pdf <- file.path(figure_dir, "fig2abc_alpha_k_benchmark_style.pdf")
out_png <- file.path(figure_dir, "fig2abc_alpha_k_benchmark_style.png")

plot_assemble(
  annotate_text("a", x = 2, y = 1, fontsize = font_size, fontface = "bold", vjust = 1),
  annotate_text("Consistency", x = 30, y = 1, fontsize = font_size_large, hjust = 0.5, vjust = 1, fontface = "bold"),
  annotate_text("b", x = 62, y = 1, fontsize = font_size, fontface = "bold", vjust = 1),
  annotate_text("Simulation", x = 90, y = 1, fontsize = font_size_large, hjust = 0.5, vjust = 1, fontface = "bold"),
  annotate_text("c", x = 122, y = 1, fontsize = font_size, fontface = "bold", vjust = 1),
  annotate_text("Downsampling", x = 150, y = 1, fontsize = font_size_large, hjust = 0.5, vjust = 1, fontface = "bold"),
  annotate_graphic(file.path(pf_dir, "illustrations", "benchmark_overviewArtboard_1.pdf"), x = 32, y = 91, units = "mm"),
  annotate_graphic(file.path(pf_dir, "illustrations", "benchmark_overviewArtboard_3.pdf"), x = 92, y = 91, units = "mm"),
  annotate_graphic(file.path(pf_dir, "illustrations", "benchmark_overviewArtboard_2.pdf"), x = 152, y = 91, units = "mm"),
  list(plot = consistency_pl, x = 0, y = 22, width = 60, height = 78),
  list(plot = simulation_pl, x = 60, y = 22, width = 60, height = 78),
  list(plot = downsampling_pl, x = 120, y = 22, width = 60, height = 78),
  width = 180,
  height = 102,
  units = "mm",
  show_grid_lines = FALSE,
  filename = out_pdf,
  latex_support = TRUE
)

if (requireNamespace("pdftools", quietly = TRUE)) {
  suppressWarnings(pdftools::pdf_convert(
    out_pdf,
    format = "png",
    dpi = 200,
    filenames = out_png
  ))
} else {
  warning("Install pdftools to also write ", out_png, call. = FALSE)
}

message("Wrote ", out_pdf)
message("Wrote ", out_png)
