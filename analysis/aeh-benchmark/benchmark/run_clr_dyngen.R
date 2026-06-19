#!/usr/bin/env Rscript
# run_clr_dyngen.R
#
# Runs CLR simulation benchmark for dyngen (5 seeds).
# Matches the original benchmark's approach:
#   - generate_dataset(config, format="sce")
#   - counts(sim) for observed UMI
#   - t(reducedDim(sim, "MDS")) as ground truth embedding
#
# Run from benchmark/ directory:
#   Rscript run_clr_dyngen.R

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(matrixStats)
  library(tibble)
  library(dyngen)
  library(SingleCellExperiment)
})

set.seed(42)

RESULTS_DIR <- "output/benchmark_results"

# ── Helpers ────────────────────────────────────────────────────────────────────
clr_transform <- function(UMI, sf) {
  if (inherits(UMI, "sparseMatrix")) UMI <- as.matrix(UMI)
  logpf <- log1p(sweep(UMI, 2L, sf, "/"))
  sweep(logpf, 2L, colMeans(logpf), "-")
}

run_pca <- function(X_genes_x_cells, n) {
  X <- t(X_genes_x_cells)
  k <- min(n, nrow(X) - 1L, ncol(X) - 1L)
  prcomp_irlba(X, n = k, center = TRUE, scale. = FALSE)$x
}

make_knn <- function(emb, k) {
  k_use <- min(k, nrow(emb) - 1L)
  FNN::get.knnx(emb, emb, k = k_use + 1L)$nn.index[, -1L, drop = FALSE][, seq_len(k_use), drop = FALSE]
}

mean_knn_overlap <- function(knn1, knn2) {
  mean(vapply(seq_len(nrow(knn1)), function(i) length(intersect(knn1[i,], knn2[i,])), 0.0))
}

size_factors <- function(UMI) { cs <- colSums(UMI); cs / mean(cs) }

filter_mat <- function(UMI, min_cell_count = 10) {
  UMI[rowSums(UMI) > 0, colSums(UMI) > min_cell_count, drop = FALSE]
}

clr_knn <- function(UMI, sf, pca_dim, knn) make_knn(run_pca(clr_transform(UMI, sf), pca_dim), knn)

SIM_KNN   <- c(10, 50, 100)
SIM_PCA   <- c(5, 10, 50)
SIM_SEEDS <- 1:5

# ── Simulate with dyngen ───────────────────────────────────────────────────────
run_dyngen_sim <- function(seed) {
  set.seed(seed)
  message("  Simulating dyngen (seed=", seed, ") ...")
  backbone <- dyngen::backbone_consecutive_bifurcating()
  num_tfs  <- nrow(backbone$module_info)
  # Match original benchmark: ~1000 total features, 5000 cells (use smaller for speed)
  # Using 500 cells, 200 features to be practical locally
  num_targets <- 100L
  num_hks     <- 50L
  config <- dyngen::initialise_model(
    backbone    = backbone,
    num_tfs     = num_tfs,
    num_targets = num_targets,
    num_hks     = num_hks,
    num_cells   = 500L,
    verbose     = FALSE
  )
  res <- dyngen::generate_dataset(config, format = "sce")
  sim <- res$dataset

  UMI  <- counts(sim)
  expressed_cells <- colSums(UMI) > 10
  expressed_genes <- rowSums(UMI) > 0
  sim  <- sim[expressed_genes, expressed_cells]
  UMI  <- counts(sim)

  # Ground truth: MDS embedding of cells (from true trajectory)
  # dim: cells x 3 → we transpose to genes_like x cells
  mds  <- reducedDim(sim, "MDS")  # cells x 3
  list(UMI = UMI, ground_truth = t(mds))  # ground_truth: 3 x ncell
}

run_sim_clr <- function(seed, sim_data) {
  if (is.null(sim_data)) return(NULL)
  UMI <- sim_data$UMI
  GT  <- sim_data$ground_truth  # 3 x ncell (MDS embedding)
  common_cells <- intersect(colnames(UMI), colnames(GT))
  if (length(common_cells) < 10) return(NULL)
  UMI <- UMI[, common_cells, drop = FALSE]
  GT  <- GT[, common_cells, drop = FALSE]
  sf  <- size_factors(UMI)

  rows <- list()
  for (pca_dim in SIM_PCA) {
    # Cache PCA per pca_dim
    pca_pred <- run_pca(clr_transform(UMI, sf), pca_dim)
    # Ground truth KNN: PCA of the MDS embedding
    gt_pca <- run_pca(as.matrix(GT), min(pca_dim, nrow(GT) - 1L))

    for (knn in SIM_KNN) {
      t0 <- proc.time()
      knn_pred <- make_knn(pca_pred, knn)
      knn_gt   <- make_knn(gt_pca, knn)
      elapsed  <- (proc.time() - t0)[["elapsed"]]
      ov <- mean_knn_overlap(knn_pred, knn_gt)
      rows[[length(rows)+1]] <- tibble(
        ARI = NA_real_, AMI = NA_real_, NMI = NA_real_,
        mean_knn_overlap = ov,
        n_clusters = NA_real_, n_clusters_counts = NA_real_,
        ground_truth_id   = paste0("gt_dyngen_s", seed),
        transformation_id = paste0("clr_dyngen_s", seed, "_p", pca_dim, "_k", knn),
        simulator = "dyngen", seed = seed, pca_dim = pca_dim, knn = knn,
        transformation = "clr", alpha = "FALSE",
        cputime_sec = elapsed, elapsed_sec = elapsed
      )
      message(sprintf("    pca=%d knn=%d  overlap=%.3f  (%.1fs)", pca_dim, knn, ov, elapsed))
    }
  }
  bind_rows(rows)
}

# ── Run ────────────────────────────────────────────────────────────────────────
message("=== DYNGEN CLR SIMULATION ===\n")
sim_rows <- list()
for (seed in SIM_SEEDS) {
  message("seed=", seed)
  sim_data <- tryCatch(run_dyngen_sim(seed),
                       error = function(e) { message("  FAILED: ", e$message); NULL })
  if (!is.null(sim_data)) {
    res <- tryCatch(run_sim_clr(seed, sim_data),
                    error = function(e) { message("  ERROR: ", e$message); NULL })
    if (!is.null(res)) sim_rows[[length(sim_rows)+1]] <- res
  }
  gc()
}

new_rows <- bind_rows(sim_rows)

if (nrow(new_rows) > 0) {
  existing <- read_tsv(file.path(RESULTS_DIR, "simulation_results.tsv"),
                       show_col_types = FALSE) %>%
    filter(!(simulator == "dyngen" & transformation == "clr")) %>%
    mutate(alpha = as.character(alpha))
  write_tsv(bind_rows(existing, new_rows),
            file.path(RESULTS_DIR, "simulation_results.tsv"))
  message("\nAdded ", nrow(new_rows), " CLR dyngen rows to simulation_results.tsv")
} else {
  message("No results computed.")
}
message("\n=== DONE ===")
