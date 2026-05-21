project_lib <- file.path(
  getwd(), "renv", "library", "R-4.5", "aarch64-apple-darwin20"
)
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(tidyverse)
library(yaml)
library(scuttle)
library(FNN)
library(MatrixGenerics)

source("src/downsampling_benchmark/download_helper.R")
source("src/transformations/transformation_helper.R")

instructions <- yaml::read_yaml("job_overview.yaml")
cfg <- instructions[["sirna_quantitative_recovery"]]

safe_cor <- function(x, y, method = "pearson") {
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = method))
}

safe_slope <- function(x, y) {
  if (length(x) < 3L || stats::sd(x) == 0) {
    return(NA_real_)
  }
  unname(stats::coef(stats::lm(y ~ x))[["x"]])
}

parse_seed_list <- function(default_seeds) {
  seed_env <- Sys.getenv("SIRNA_SEEDS", unset = "")
  if (seed_env == "") {
    return(default_seeds)
  }

  as.integer(trimws(strsplit(seed_env, ",", fixed = TRUE)[[1]]))
}

signature_from_training <- function(expr_mat, group_idx, ctrl_idx, n_genes) {
  mean_group <- MatrixGenerics::rowMeans2(expr_mat[, group_idx, drop = FALSE])
  mean_ctrl <- MatrixGenerics::rowMeans2(expr_mat[, ctrl_idx, drop = FALSE])
  effect <- mean_group - mean_ctrl
  effect[!is.finite(effect)] <- 0

  list(
    up = head(order(effect, decreasing = TRUE), n_genes),
    down = head(order(effect, decreasing = FALSE), n_genes)
  )
}

score_signature <- function(expr_mat, signature) {
  MatrixGenerics::colMeans2(expr_mat[signature$up, , drop = FALSE]) -
    MatrixGenerics::colMeans2(expr_mat[signature$down, , drop = FALSE])
}

compute_seed_dataset <- function(seed) {
  set.seed(seed)
  sce <- data_loaders[["smartSeq3_siRNA_knockdown"]]()
  UMI <- as.matrix(SummarizedExperiment::assay(sce))
  col_data <- as.data.frame(SummarizedExperiment::colData(sce))

  colsums <- colSums2(UMI)
  downsample_proportion <- 5000 / stats::median(colsums)
  reduced <- as.matrix(scuttle::downsampleMatrix(UMI, prop = downsample_proportion, bycol = FALSE))

  expressed_cells <- matrixStats::colSums2(reduced) > 0
  expressed_genes <- matrixStats::rowSums2(reduced) > 0

  list(
    full = UMI[expressed_genes, expressed_cells, drop = FALSE],
    reduced = reduced[expressed_genes, expressed_cells, drop = FALSE],
    col_data = col_data[expressed_cells, , drop = FALSE]
  )
}

run_one_seed <- function(seed, pca_dim, regression_knn, trans_specs, signature_n_genes = 50L) {
  dataset <- compute_seed_dataset(seed)
  labels <- as.character(dataset$col_data$predicted)

  sf_full <- MatrixGenerics::colSums2(dataset$full)
  sf_full <- sf_full / mean(sf_full)
  reference_expr <- logp1_fnc(dataset$full, sf_full, FALSE)

  reps <- list()
  for (spec in trans_specs) {
    trans_name <- spec$name
    alpha <- spec$alpha[[1]]
    alpha_parsed <- if (alpha == "TRUE") TRUE else if (alpha == "FALSE") FALSE else readr::parse_double(alpha)
    sf_reduced <- MatrixGenerics::colSums2(dataset$reduced)
    sf_reduced <- sf_reduced / mean(sf_reduced)
    message(sprintf("seed=%d transform=%s", seed, trans_name))
    trans_dat <- all_transformations[[trans_name]](dataset$reduced, sf_reduced, alpha_parsed)
    reps[[trans_name]] <- list(
      rep = make_lowdim_representation(trans_name, trans_dat, pca_dim),
      alpha = alpha
    )
  }

  control_label <- "Cont"
  perturbations <- sort(setdiff(unique(labels), control_label))
  results <- list()

  set.seed(seed)
  for (perturbation in perturbations) {
    group_idx <- which(labels == perturbation)
    ctrl_idx <- which(labels == control_label)
    train_group <- sample(group_idx, floor(length(group_idx) * 0.7))
    train_ctrl <- sample(ctrl_idx, floor(length(ctrl_idx) * 0.7))
    test_group <- setdiff(group_idx, train_group)

    signature <- signature_from_training(reference_expr, train_group, train_ctrl, signature_n_genes)
    target_score <- score_signature(reference_expr, signature)
    train_idx <- c(train_group, train_ctrl)

    for (trans_name in names(reps)) {
      rep <- reps[[trans_name]]$rep
      k_use <- min(regression_knn, length(train_idx))
      pred <- FNN::knn.reg(
        train = rep[train_idx, , drop = FALSE],
        test = rep[test_group, , drop = FALSE],
        y = target_score[train_idx],
        k = k_use
      )$pred
      truth <- target_score[test_group]

      results[[length(results) + 1L]] <- tibble(
        dataset = "smartSeq3_siRNA_knockdown",
        seed = seed,
        pca_dim = pca_dim,
        perturbation = perturbation,
        n_train = length(train_idx),
        n_test = length(test_group),
        regression_knn = k_use,
        signature_n_genes = signature_n_genes,
        transformation = trans_name,
        alpha = reps[[trans_name]]$alpha,
        pearson = safe_cor(truth, pred, method = "pearson"),
        spearman = safe_cor(truth, pred, method = "spearman"),
        rmse = sqrt(mean((pred - truth) ^ 2)),
        slope = safe_slope(truth, pred),
        sd_ratio = if (stats::sd(truth) == 0) NA_real_ else stats::sd(pred) / stats::sd(truth)
      )
    }
  }

  bind_rows(results)
}

dir.create("output/benchmark_results", recursive = TRUE, showWarnings = FALSE)

requested_transforms <- c(
  "logp1",
  "logp1_hvg",
  "logp1_zscore",
  "logp1_hvg_zscore",
  "logp_cpm",
  "logp1_size_normed",
  "acosh",
  "logp_alpha",
  "pearson",
  "pearson_clip",
  "pearson_clip_hvg",
  "pearson_clip_zscore",
  "pearson_clip_hvg_zscore",
  "scgpt",
  "clr"
)

trans_specs <- cfg$representation_construction$transformations
trans_specs <- purrr::keep(trans_specs, ~ .x$name %in% requested_transforms)
seeds <- parse_seed_list(cfg$input_data$seed)
pca_dim <- cfg$representation_construction$pca[[1]]
regression_knn <- cfg$representation_construction$regression_knn[[1]]
result_file <- Sys.getenv(
  "SIRNA_RESULT_FILE",
  unset = "output/benchmark_results/sirna_quantitative_recovery_results.tsv"
)

if (file.exists(result_file)) {
  all_results <- readr::read_tsv(
    result_file,
    show_col_types = FALSE,
    col_types = readr::cols(alpha = readr::col_character())
  ) %>%
    mutate(alpha = as.character(alpha))
} else {
  all_results <- tibble()
}

for (seed in seeds) {
  if ("seed" %in% colnames(all_results) && seed %in% all_results$seed) {
    message(sprintf("Skipping seed=%d because it already exists in %s", seed, result_file))
    next
  }

  seed_results <- run_one_seed(seed, pca_dim, regression_knn, trans_specs)
  seed_results <- seed_results %>%
    mutate(alpha = as.character(alpha))
  all_results <- bind_rows(all_results, seed_results) %>%
    distinct(seed, perturbation, transformation, alpha, .keep_all = TRUE) %>%
    arrange(seed, transformation, perturbation)

  readr::write_tsv(all_results, result_file)
  message(sprintf("Saved through seed=%d to %s", seed, result_file))
}

message(sprintf("Finished %d seed(s); results in %s", length(seeds), result_file))
