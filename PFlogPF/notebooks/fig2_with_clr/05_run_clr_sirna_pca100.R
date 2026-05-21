#!/usr/bin/env Rscript
# 05_run_clr_sirna_pca100.R
#
# Computes CLR downsampling benchmark for smartSeq3_siRNA_knockdown at
# pca_dim = 100 (extends the default set of 5, 10, 50 for the PCA-dim plot).
#
# Run from notebooks/fig2_with_clr/ directory:
#   Rscript 05_run_clr_sirna_pca100.R

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(tibble)
})

RESULTS_DIR <- "../../benchmark/output/benchmark_results"
DATA_DIR    <- "../../benchmark/output/clr_local/data/downsampling"

clr_transform <- function(UMI, sf) {
  if (inherits(UMI, "sparseMatrix")) UMI <- as.matrix(UMI)
  logpf <- log1p(sweep(UMI, 2L, sf, "/"))
  sweep(logpf, 2L, colMeans(logpf), "-")
}
run_pca <- function(X, n) {
  X <- t(X); k <- min(n, nrow(X) - 1L, ncol(X) - 1L)
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
filter_matrix <- function(M) M[rowSums(M) > 0, colSums(M) > 0, drop = FALSE]

# ── Load siRNA data ───────────────────────────────────────────────────────────
f_cast <- file.path(DATA_DIR, "ss3_n4298_fibs_siKD_umiCast.rds")
f_c57  <- file.path(DATA_DIR, "ss3_n4298_fibs_siKD_umiC57.rds")
cast_raw <- readRDS(f_cast); c57_raw <- readRDS(f_c57)
total <- cast_raw + c57_raw; total[is.na(total)] <- 0L
UMI_full <- filter_matrix(as.matrix(total))
message(ncol(UMI_full), " cells x ", nrow(UMI_full), " genes")

# ── Run pca_dim = 100 for seeds 1:5, knn 10/50/100 ───────────────────────────
pca_dim <- 100L
rows <- list()
for (seed in 1:5) {
  set.seed(seed)
  UMI_red <- apply(UMI_full, 2, function(cell) {
    n <- max(1L, round(sum(cell) * 0.1))
    as.integer(rmultinom(1L, size = n, prob = cell + 1e-8))
  })
  rownames(UMI_red) <- rownames(UMI_full)
  UMI_red <- filter_matrix(UMI_red)
  common  <- intersect(rownames(UMI_full), rownames(UMI_red))
  UMI_f <- UMI_full[common, , drop = FALSE]; UMI_r <- UMI_red[common, , drop = FALSE]
  sf_f  <- size_factors(UMI_f);              sf_r  <- size_factors(UMI_r)
  pca_f <- run_pca(clr_transform(UMI_f, sf_f), pca_dim)
  pca_r <- run_pca(clr_transform(UMI_r, sf_r), pca_dim)
  for (knn in c(10, 50, 100)) {
    t0 <- proc.time()
    knn_f <- make_knn(pca_f, knn); knn_r <- make_knn(pca_r, knn)
    elapsed <- (proc.time() - t0)[["elapsed"]]
    ov <- mean_knn_overlap(knn_f, knn_r)
    rows[[length(rows) + 1]] <- tibble(
      overlap = ov,
      transformation_full_data_ids    = paste0("clr_local_full_smartSeq3_siRNA_knockdown_s", seed, "_p", pca_dim, "_k", knn),
      transformation_reduced_data_ids = paste0("clr_local_red_smartSeq3_siRNA_knockdown_s",  seed, "_p", pca_dim, "_k", knn),
      dataset = "smartSeq3_siRNA_knockdown", seed = seed, pca_dim = pca_dim, knn = knn,
      transformation = "clr", alpha = "FALSE",
      full_cputime_sec = elapsed / 2, full_elapsed_sec = elapsed / 2,
      reduced_cputime_sec = elapsed / 2, reduced_elapsed_sec = elapsed / 2
    )
    message(sprintf("seed=%d knn=%d  overlap=%.4f  (%.1fs)", seed, knn, ov, elapsed))
  }
}

new_rows <- bind_rows(rows)
existing <- read_tsv(file.path(RESULTS_DIR, "downsampling_results.tsv"),
                     show_col_types = FALSE) %>%
  mutate(alpha = as.character(alpha)) %>%
  filter(!(dataset == "smartSeq3_siRNA_knockdown" & transformation == "clr" & pca_dim == 100))
write_tsv(bind_rows(existing, new_rows),
          file.path(RESULTS_DIR, "downsampling_results.tsv"))
message("Added ", nrow(new_rows), " rows to downsampling_results.tsv")
message("=== DONE ===")
