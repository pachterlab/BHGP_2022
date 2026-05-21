#!/usr/bin/env Rscript
# run_clr_downsampling_missing.R
#
# Runs CLR downsampling benchmark for the 2 missing datasets:
#   - smartSeq3_hek  (Smartseq3.HEK.cleanup.UMIcounts.txt from E-MTAB-8735)
#   - smartSeq3_siRNA_knockdown (RDS files from sandberg-lab GitHub)
#
# Run from benchmark/ directory:
#   Rscript run_clr_downsampling_missing.R

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(matrixStats)
  library(tibble)
  library(purrr)
})

set.seed(42)

RESULTS_DIR <- "output/benchmark_results"
DATA_DIR    <- "output/clr_local/data/downsampling"

# ── Shared helpers ─────────────────────────────────────────────────────────────
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
  mean(vapply(seq_len(nrow(knn1)), function(i)
    length(intersect(knn1[i,], knn2[i,])), 0.0))
}

size_factors <- function(UMI) { cs <- colSums(UMI); cs / mean(cs) }

filter_matrix <- function(M) M[rowSums(M) > 0, colSums(M) > 0, drop = FALSE]

clr_knn <- function(UMI, sf, pca_dim, knn) make_knn(run_pca(clr_transform(UMI, sf), pca_dim), knn)

DS_KNN   <- c(10, 50, 100)
DS_PCA   <- c(5, 10, 50)
DS_SEEDS <- 1:5

# ── Data loaders ───────────────────────────────────────────────────────────────
load_hek <- function() {
  f <- file.path(DATA_DIR, "Smartseq3.HEK.cleanup.UMIcounts.txt")
  stopifnot(file.exists(f))
  as.matrix(read.delim(f))
}

load_sirna <- function() {
  f_cast <- file.path(DATA_DIR, "ss3_n4298_fibs_siKD_umiCast.rds")
  f_c57  <- file.path(DATA_DIR, "ss3_n4298_fibs_siKD_umiC57.rds")
  stopifnot(file.exists(f_cast), file.exists(f_c57))
  cast_raw <- readRDS(f_cast)
  c57_raw  <- readRDS(f_c57)
  stopifnot(all(rownames(cast_raw) == rownames(c57_raw)),
            all(colnames(cast_raw) == colnames(c57_raw)))
  total <- cast_raw + c57_raw
  # Replace NAs with 0 (imputed positions)
  total[is.na(total)] <- 0L
  as.matrix(total)
}

# ── Run downsampling CLR ───────────────────────────────────────────────────────
run_downsampling_clr <- function(dataset_name, UMI_full) {
  UMI_full <- filter_matrix(UMI_full)
  message("  ", dataset_name, ": ", ncol(UMI_full), " cells x ", nrow(UMI_full), " genes")

  rows <- list()
  for (seed in DS_SEEDS) {
    set.seed(seed)
    UMI_red <- apply(UMI_full, 2, function(cell) {
      n <- max(1L, round(sum(cell) * 0.1))
      as.integer(rmultinom(1L, size = n, prob = cell + 1e-8))
    })
    rownames(UMI_red) <- rownames(UMI_full)
    UMI_red <- filter_matrix(UMI_red)
    common_genes <- intersect(rownames(UMI_full), rownames(UMI_red))
    UMI_f  <- UMI_full[common_genes, , drop = FALSE]
    UMI_r  <- UMI_red[common_genes, , drop = FALSE]
    sf_f   <- size_factors(UMI_f)
    sf_r   <- size_factors(UMI_r)

    for (pca_dim in DS_PCA) {
      for (knn in DS_KNN) {
        t0 <- proc.time()
        knn_full <- clr_knn(UMI_f, sf_f, pca_dim, knn)
        knn_red  <- clr_knn(UMI_r, sf_r, pca_dim, knn)
        elapsed  <- (proc.time() - t0)[["elapsed"]]
        ov <- mean_knn_overlap(knn_full, knn_red)
        rows[[length(rows)+1]] <- tibble(
          overlap = ov,
          transformation_full_data_ids    = paste0("clr_local_full_", dataset_name, "_s", seed, "_p", pca_dim, "_k", knn),
          transformation_reduced_data_ids = paste0("clr_local_red_", dataset_name, "_s", seed, "_p", pca_dim, "_k", knn),
          dataset = dataset_name, seed = seed, pca_dim = pca_dim, knn = knn,
          transformation = "clr", alpha = "FALSE",
          full_cputime_sec = elapsed / 2, full_elapsed_sec = elapsed / 2,
          reduced_cputime_sec = elapsed / 2, reduced_elapsed_sec = elapsed / 2
        )
        message(sprintf("    seed=%d pca=%d knn=%d  overlap=%.4f  (%.1fs)",
                        seed, pca_dim, knn, ov, elapsed))
      }
    }
  }
  bind_rows(rows)
}

# ── Main ───────────────────────────────────────────────────────────────────────
message("=== DOWNSAMPLING CLR — MISSING DATASETS ===\n")

results <- list()

message("Dataset: smartSeq3_hek")
hek_mat <- tryCatch(load_hek(), error = function(e) { message("FAILED: ", e$message); NULL })
if (!is.null(hek_mat)) {
  res <- tryCatch(run_downsampling_clr("smartSeq3_hek", hek_mat),
                  error = function(e) { message("ERROR: ", e$message); NULL })
  if (!is.null(res)) results[[1]] <- res
}

message("\nDataset: smartSeq3_siRNA_knockdown")
sirna_mat <- tryCatch(load_sirna(), error = function(e) { message("FAILED: ", e$message); NULL })
if (!is.null(sirna_mat)) {
  res <- tryCatch(run_downsampling_clr("smartSeq3_siRNA_knockdown", sirna_mat),
                  error = function(e) { message("ERROR: ", e$message); NULL })
  if (!is.null(res)) results[[2]] <- res
}

new_rows <- bind_rows(results)

if (nrow(new_rows) > 0) {
  existing <- read_tsv(file.path(RESULTS_DIR, "downsampling_results.tsv"),
                       show_col_types = FALSE) %>%
    mutate(alpha = as.character(alpha))
  write_tsv(bind_rows(existing, new_rows),
            file.path(RESULTS_DIR, "downsampling_results.tsv"))
  message("\nAdded ", nrow(new_rows), " CLR rows to downsampling_results.tsv")
} else {
  message("No new results computed.")
}

message("\n=== DONE ===")
