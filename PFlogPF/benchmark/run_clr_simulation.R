#!/usr/bin/env Rscript
# run_clr_simulation.R
# Runs CLR and CLR_alpha through the simulation benchmark and appends results
# to simulation_results.tsv.
#
# Run from the benchmark/ directory:
#   Rscript run_clr_simulation.R

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(SingleCellExperiment)
  library(BiocNeighbors)
  library(bluster)
  library(igraph)
  library(aricode)
})

RESULTS_DIR <- file.path(getwd(), "output/benchmark_results")
DATA_DIR    <- file.path(getwd(), "output/clr_local/data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

SEEDS    <- 1:5
KNN_VALS <- c(10, 50, 100)
PCA_VALS <- c(5, 10, 50)
# Skip dyngen (very slow). muscat, linear_walk, random_walk, scDesign2.
SIMULATORS <- c("muscat", "linear_walk", "random_walk", "scDesign2")

# ── Check what's already done ─────────────────────────────────────────────────
sim_tsv <- file.path(RESULTS_DIR, "simulation_results.tsv")
existing <- read_tsv(sim_tsv, show_col_types = FALSE) %>%
  filter(transformation %in% c("clr", "clr_alpha"))
cat("Existing CLR/CLR_alpha simulation rows:", nrow(existing), "\n")

# ── Transforms ────────────────────────────────────────────────────────────────
clr_transform <- function(UMI) {
  sf <- Matrix::colSums(UMI); sf <- sf / mean(sf)
  logpf <- log1p(sweep(as.matrix(UMI), 2L, sf, "/"))
  sweep(logpf, 2L, colMeans(logpf), "-")
}
clr_alpha_transform <- function(UMI, alpha = 0.05) {
  sf <- Matrix::colSums(UMI); sf <- sf / mean(sf)
  logpf <- log(sweep(as.matrix(UMI), 2L, sf, "/") + 1 / (4 * alpha))
  sweep(logpf, 2L, colMeans(logpf), "-")
}

run_pca <- function(mat_genes_x_cells, n) {
  X <- t(mat_genes_x_cells)
  k <- min(n, nrow(X) - 1L, ncol(X) - 1L)
  prcomp_irlba(X, n = k, center = TRUE, scale. = FALSE)$x
}

make_knn <- function(pca, k) {
  nn <- FNN::get.knnx(pca, pca, k = k + 1L)$nn.index
  nn[, -1L, drop = FALSE]
}

# Ground-truth KNN using Annoy (same as original benchmark)
ground_truth_knn <- function(gt_mat, k) {
  # gt_mat is genes x cells; transpose for BiocNeighbors
  BiocNeighbors::findAnnoy(t(gt_mat), k = k, warn.ties = FALSE,
                            get.distance = FALSE)$index
}

compute_overlap <- function(knn_pred, knn_gt) {
  mean(vapply(seq_len(nrow(knn_pred)), function(i)
    length(intersect(knn_pred[i,], knn_gt[i,])), 0.0))
}

compute_clustering_metrics <- function(knn_pred, knn_gt, k) {
  gt_graph  <- bluster::neighborsToKNNGraph(knn_gt)
  gt_clust  <- factor(igraph::membership(igraph::cluster_walktrap(gt_graph)))
  pred_graph <- bluster::neighborsToKNNGraph(knn_pred)
  pred_clust <- factor(igraph::membership(igraph::cluster_walktrap(pred_graph)))
  list(
    ARI = aricode::ARI(gt_clust, pred_clust),
    AMI = aricode::AMI(gt_clust, pred_clust),
    NMI = aricode::NMI(gt_clust, pred_clust),
    n_clusters = nlevels(gt_clust),
    n_clusters_counts = nlevels(pred_clust)
  )
}

# ── Simulators ────────────────────────────────────────────────────────────────

simulate_muscat <- function(seed) {
  set.seed(seed)
  library(muscat)
  data(example_sce)
  sce_preped <- prepSim(example_sce, verbose = FALSE)
  sim <- muscat::simData(sce_preped, rel_lfc = c(1, 0.5, 0.1, 0.05),
                         nc = 5e3, nk = 4,
                         p_dd = c(0.7, 0, 0.3, 0, 0, 0),
                         lfc = 2, ng = 1e3, force = TRUE)
  tmp <- as_tibble(metadata(sim)$gene_info) %>%
    pivot_longer(starts_with("sim_mean"), names_sep = "\\.",
                 names_to = c(".value", "group_id")) %>%
    dplyr::select(gene, cluster_id, group_id, sim_mean) %>%
    full_join(colData(sim) %>% as_tibble(rownames = "cell_id") %>%
                mutate(cell_id = factor(cell_id, levels = cell_id)),
              by = c("cluster_id", "group_id"))
  gt <- tmp %>% arrange(cell_id) %>%
    dplyr::select(gene, cell_id, sim_mean) %>%
    pivot_wider(id_cols = gene, names_from = cell_id, values_from = sim_mean) %>%
    column_to_rownames("gene") %>% as.matrix()
  list(ground_truth = log10(gt + 1e-4), UMI = assay(sim, "counts"))
}

simulate_linear_walk <- function(seed) {
  set.seed(seed)
  sce <- scRNAseq::BaronPancreasData("human")
  sce <- scuttle::logNormCounts(sce)
  ref <- as.matrix(assay(sce, "counts"))
  ref <- ref[Matrix::rowSums(ref) > 0, Matrix::colSums(ref) > 0]
  n_genes <- nrow(ref); n_cells <- ncol(ref)
  N_c <- colSums(ref)
  mu_tilde_g <- log(rowSums(ref) / sum(ref))
  sig2_g <- rexp(n_genes, rate = 1/2)
  branch_length <- 1200
  delta_true <- matrix(NA, n_genes, n_cells)
  branch_idx <- 0; start_point <- NULL; end_point <- NULL
  for (idx in seq_len(n_cells) - 1) {
    if (idx == 0) {
      start_point <- rnorm(n_genes); end_point <- rnorm(n_genes)
      delta_true[, 1] <- rnorm(n_genes, start_point, 0.1)
    } else if (idx %% branch_length == 0) {
      branch_idx <- 0
      start_point <- rnorm(n_genes, delta_true[, sample.int(idx,1)], 0.1)
      end_point   <- rnorm(n_genes, 0, 5)
      delta_true[, idx+1] <- rnorm(n_genes,
        start_point + (end_point - start_point) * branch_idx / branch_length, 0.1)
    } else {
      branch_idx <- branch_idx + 1
      delta_true[, idx+1] <- rnorm(n_genes,
        start_point + (end_point - start_point) * branch_idx / branch_length, 0.1)
    }
  }
  lambda <- sqrt(sig2_g / rowVars(delta_true))
  delta_true <- (delta_true - rowMeans(delta_true)) * lambda
  mu_g <- mu_tilde_g - sig2_g / 2
  UMI <- matrix(rnbinom(n_genes * n_cells,
                         mu = t(t(exp(mu_g + delta_true)) * N_c), size = 100),
                n_genes, n_cells)
  list(ground_truth = delta_true, UMI = UMI)
}

simulate_random_walk <- function(seed) {
  set.seed(seed)
  sce <- scRNAseq::BaronPancreasData("human")
  sce <- scuttle::logNormCounts(sce)
  ref <- as.matrix(assay(sce, "counts"))
  ref <- ref[Matrix::rowSums(ref) > 0, Matrix::colSums(ref) > 0]
  n_genes <- nrow(ref); n_cells <- ncol(ref)
  N_c <- colSums(ref)
  mu_tilde_g <- log(rowSums(ref) / sum(ref))
  sig2_g <- rexp(n_genes, rate = 1/2)
  branch_length <- 13
  delta_true <- matrix(NA, n_genes, n_cells)
  parents <- rep(NA, n_cells)
  for (idx in seq_len(n_cells)) {
    if (idx == 1) {
      delta_true[, 1] <- rnorm(n_genes); parents[1] <- 0
    } else if (idx %% branch_length == 0) {
      parent_id <- sample.int(idx - 1, 1)
      parents[idx] <- parent_id
      delta_true[, idx] <- rnorm(n_genes, delta_true[, parent_id], 1)
    } else {
      parents[idx] <- idx - 1
      delta_true[, idx] <- rnorm(n_genes, delta_true[, idx - 1], 1)
    }
  }
  lambda <- sqrt(sig2_g / rowVars(delta_true))
  delta_true <- (delta_true - rowMeans(delta_true)) * lambda
  mu_g <- mu_tilde_g - sig2_g / 2
  UMI <- matrix(rnbinom(n_genes * n_cells,
                         mu = t(t(exp(mu_g + delta_true)) * N_c), size = 100),
                n_genes, n_cells)
  list(ground_truth = delta_true, UMI = UMI)
}

simulate_scDesign2 <- function(seed) {
  set.seed(seed)
  # Download GSE130931 if not cached
  gse_ids <- c("GSM4041124", "GSM4041125")
  gse_dir <- file.path(DATA_DIR, "consistency")
  dirs <- file.path(gse_dir, gse_ids)
  if (!all(file.exists(file.path(dirs, "matrix.mtx.gz")))) {
    for (id in gse_ids) {
      out <- file.path(gse_dir, id)
      dir.create(out, recursive = TRUE, showWarnings = FALSE)
      GEOquery::getGEOSuppFiles(id, baseDir = gse_dir, makeDirectory = TRUE)
      for (f in list.files(out)) {
        nf <- f
        nf <- sub("^.*matrix\\.mtx",   "matrix.mtx",   nf)
        nf <- sub("^.*barcodes\\.tsv", "barcodes.tsv", nf)
        nf <- sub("^.*features\\.tsv", "features.tsv", nf)
        nf <- sub("^.*genes\\.tsv",    "genes.tsv",    nf)
        if (nf != f) file.rename(file.path(out, f), file.path(out, nf))
      }
    }
  }
  sce <- DropletUtils::read10xCounts(dirs)
  colData(sce)$cluster_id <- scran::quickCluster(sce, min.size = 20)
  mat <- as.matrix(assay(sce, "counts"))
  mat <- mat[Matrix::rowSums(mat) > 0, Matrix::colSums(mat) > 0]
  cls <- as.character(colData(sce)$cluster_id[Matrix::colSums(assay(sce)) > 0])
  colnames(mat) <- cls
  cluster_names <- unique(cls)
  fit <- scDesign2::fit_model_scDesign2(mat, cell_type_sel = cluster_names,
                                         sim_method = "copula", marginal = "nb")
  UMI <- scDesign2::simulate_count_scDesign2(
    fit, n_cell_new = ncol(mat),
    cell_type_prop = table(cls)[cluster_names] / length(cls),
    sim_method = "copula")
  gt <- do.call(cbind, lapply(fit[cluster_names], function(fi) {
    m <- matrix(NA, nrow = nrow(mat), ncol = fi$n_cell)
    m[fi$gene_sel1, ] <- fi$marginal_param1[, 3]
    m[fi$gene_sel2, ] <- fi$marginal_param2[, 3]
    m[fi$gene_sel3, ] <- 1e-8
    m
  }))
  list(ground_truth = log10(gt + 1e-8), UMI = UMI)
}

SIMULATOR_FNS <- list(
  muscat      = simulate_muscat,
  linear_walk = simulate_linear_walk,
  random_walk = simulate_random_walk,
  scDesign2   = simulate_scDesign2
)

# ── Main loop ─────────────────────────────────────────────────────────────────
results <- list()

for (sim_name in SIMULATORS) {
  cat("\n══ Simulator:", sim_name, "══\n")
  sim_fn <- SIMULATOR_FNS[[sim_name]]

  for (seed in SEEDS) {
    # Check if already done
    done <- existing %>%
      filter(simulator == sim_name, seed == !!seed)
    if (nrow(done) == length(KNN_VALS) * length(PCA_VALS) * 2) {  # 2 transforms
      cat("  seed", seed, "— already complete, skipping\n")
      next
    }

    cat("  seed", seed, "— simulating... ")
    sim_data <- tryCatch(sim_fn(seed), error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n"); NULL
    })
    if (is.null(sim_data)) next
    cat("done (", nrow(sim_data$UMI), "genes x", ncol(sim_data$UMI), "cells)\n")

    UMI <- sim_data$UMI
    gt  <- sim_data$ground_truth
    # Filter zeros
    keep_g <- Matrix::rowSums(UMI) > 0
    keep_c <- Matrix::colSums(UMI) > 0
    UMI <- UMI[keep_g, keep_c]
    gt  <- gt[keep_g, keep_c]

    # Apply transforms
    cat("    transforming... ")
    trans_list <- list(
      clr       = tryCatch(clr_transform(UMI),       error=function(e) NULL),
      clr_alpha = tryCatch(clr_alpha_transform(UMI), error=function(e) NULL)
    )
    cat("done\n")

    for (pca_dim in PCA_VALS) {
      cat("    pca =", pca_dim, "... ")
      for (trans_name in names(trans_list)) {
        trans <- trans_list[[trans_name]]
        if (is.null(trans)) next

        pca <- tryCatch(run_pca(trans, pca_dim), error = function(e) NULL)
        if (is.null(pca)) next

        for (k in KNN_VALS) {
          knn_pred <- tryCatch(make_knn(pca, k), error = function(e) NULL)
          if (is.null(knn_pred)) next

          knn_gt <- tryCatch(ground_truth_knn(gt, k), error = function(e) NULL)
          if (is.null(knn_gt)) next

          ov  <- compute_overlap(knn_pred, knn_gt)
          cls <- tryCatch(compute_clustering_metrics(knn_pred, knn_gt, k),
                          error = function(e) list(ARI=NA,AMI=NA,NMI=NA,
                                                   n_clusters=NA,n_clusters_counts=NA))

          results[[length(results) + 1]] <- tibble(
            ARI                = cls$ARI,
            AMI                = cls$AMI,
            NMI                = cls$NMI,
            mean_knn_overlap   = ov,
            n_clusters         = cls$n_clusters,
            n_clusters_counts  = cls$n_clusters_counts,
            ground_truth_id    = "local",
            transformation_id  = "local",
            simulator          = sim_name,
            seed               = seed,
            pca_dim            = pca_dim,
            knn                = k,
            transformation     = trans_name,
            alpha              = "FALSE",
            cputime_sec        = NA_real_,
            elapsed_sec        = NA_real_
          )
        }
      }
      cat("done\n")
    }
  }
}

# ── Append results ─────────────────────────────────────────────────────────────
if (length(results) > 0) {
  new_rows <- bind_rows(results)
  cat("\nNew rows computed:", nrow(new_rows), "\n")

  existing_tsv <- read_tsv(sim_tsv, show_col_types = FALSE)
  # Remove any prior clr rows (avoid duplicates on re-run)
  existing_tsv <- existing_tsv %>%
    filter(!(transformation %in% c("clr","clr_alpha") & simulator %in% SIMULATORS))
  combined <- bind_rows(existing_tsv, new_rows)
  write_tsv(combined, sim_tsv)
  cat("Appended to", sim_tsv, "\n")
  cat("\nSummary by simulator/transformation:\n")
  print(new_rows %>% group_by(simulator, transformation) %>%
          summarize(n=n(), mean_overlap=round(mean(mean_knn_overlap,na.rm=TRUE),3),
                    mean_ARI=round(mean(ARI,na.rm=TRUE),3), .groups="drop"))
} else {
  cat("No new results computed.\n")
}
