library(tidyverse)
library(FNN)

pa <- argparser::arg_parser("Compute quantitative perturbation-recovery metrics for the siRNA benchmark")
pa <- argparser::add_argument(pa, "--representation_result_ids", type = "character", nargs = Inf, help = "A list of reduced-data representation result ids")
pa <- argparser::add_argument(pa, "--input_dataset_id", type = "character", help = "The dataset object id in output/results")
pa <- argparser::add_argument(pa, "--dataset", type = "character", help = "[Just for documentation purposes] A readable identifier of the data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "[Just for documentation purposes] The seed used to tame randomness")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "[Just for documentation purposes] The number of dimensions in the representation")
pa <- argparser::add_argument(pa, "--transformations", type = "character", nargs = Inf, help = "[Just for documentation purposes] A readable identifier of the transformation")
pa <- argparser::add_argument(pa, "--alphas", type = "character", nargs = Inf, help = "[Just for documentation purposes] Specification of the overdispersion.")
pa <- argparser::add_argument(pa, "--regression_knn", type = "numeric", default = 30, help = "The number of neighbors used in kNN regression")
pa <- argparser::add_argument(pa, "--signature_n_genes", type = "numeric", default = 50, help = "Top genes per direction used in the perturbation signature")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)
stopifnot(pa$dataset == "smartSeq3_siRNA_knockdown")
stopifnot(length(pa$transformations) == length(pa$representation_result_ids))
stopifnot(length(pa$alphas) == length(pa$representation_result_ids))

source("src/transformations/transformation_helper.R")

safe_cor <- function(x, y, method = "pearson"){
  if(length(x) < 3L || sd(x) == 0 || sd(y) == 0){
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, method = method))
}

safe_slope <- function(x, y){
  if(length(x) < 3L || sd(x) == 0){
    return(NA_real_)
  }
  coef(lm(y ~ x))[["x"]]
}

signature_from_training <- function(expr_mat, group_idx, ctrl_idx, n_genes){
  mean_group <- MatrixGenerics::rowMeans2(expr_mat[, group_idx, drop = FALSE])
  mean_ctrl <- MatrixGenerics::rowMeans2(expr_mat[, ctrl_idx, drop = FALSE])
  effect <- mean_group - mean_ctrl
  effect[!is.finite(effect)] <- 0

  up <- head(order(effect, decreasing = TRUE), n_genes)
  down <- head(order(effect, decreasing = FALSE), n_genes)

  list(up = up, down = down)
}

score_signature <- function(expr_mat, signature){
  up_score <- MatrixGenerics::colMeans2(expr_mat[signature$up, , drop = FALSE])
  down_score <- MatrixGenerics::colMeans2(expr_mat[signature$down, , drop = FALSE])
  up_score - down_score
}

dataset <- readRDS(file.path(pa$working_dir, "results", pa$input_dataset_id))
labels <- as.character(dataset$col_data$predicted)
stopifnot(length(labels) == ncol(dataset$full))
stopifnot(identical(labels, as.character(dataset$col_data$predicted)))

full_mat <- dataset$full
sf_full <- MatrixGenerics::colSums2(full_mat)
sf_full <- sf_full / mean(sf_full)
reference_expr <- logp1_fnc(full_mat, sf_full, FALSE)

reps <- lapply(pa$representation_result_ids, function(id){
  rep <- readRDS(file.path(pa$working_dir, "results", id))
  stopifnot(nrow(rep) == ncol(full_mat))
  as.matrix(rep)
})

control_label <- "Cont"
perturbations <- sort(setdiff(unique(labels), control_label))

set.seed(pa$seed)
results <- list()

for(perturbation in perturbations){
  group_idx <- which(labels == perturbation)
  ctrl_idx <- which(labels == control_label)
  train_group <- sample(group_idx, floor(length(group_idx) * 0.7))
  train_ctrl <- sample(ctrl_idx, floor(length(ctrl_idx) * 0.7))
  test_group <- setdiff(group_idx, train_group)

  signature <- signature_from_training(reference_expr, train_group, train_ctrl, pa$signature_n_genes)
  target_score <- score_signature(reference_expr, signature)
  train_idx <- c(train_group, train_ctrl)

  for(i in seq_along(reps)){
    rep <- reps[[i]]
    k_use <- min(pa$regression_knn, length(train_idx))
    pred <- FNN::knn.reg(
      train = rep[train_idx, , drop = FALSE],
      test = rep[test_group, , drop = FALSE],
      y = target_score[train_idx],
      k = k_use
    )$pred

    truth <- target_score[test_group]
    results[[length(results) + 1L]] <- tibble(
      dataset = pa$dataset,
      seed = pa$seed,
      pca_dim = pa$pca_dim,
      perturbation = perturbation,
      n_train = length(train_idx),
      n_test = length(test_group),
      regression_knn = k_use,
      signature_n_genes = pa$signature_n_genes,
      transformation = pa$transformations[[i]],
      alpha = pa$alphas[[i]],
      representation_result_id = pa$representation_result_ids[[i]],
      pearson = safe_cor(truth, pred, method = "pearson"),
      spearman = safe_cor(truth, pred, method = "spearman"),
      rmse = sqrt(mean((pred - truth) ^ 2)),
      slope = safe_slope(truth, pred),
      sd_ratio = if(sd(truth) == 0) NA_real_ else sd(pred) / sd(truth)
    )
  }
}

write_tsv(bind_rows(results), file.path(pa$working_dir, "results", pa$result_id))
