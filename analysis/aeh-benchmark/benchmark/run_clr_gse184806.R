#!/usr/bin/env Rscript
# run_clr_gse184806.R
#
# Runs CLR consistency benchmark for GSE184806 (43k cells) using a
# memory-efficient approach that avoids forming the full dense CLR matrix.
#
# The centered CLR matrix is: ctA[j,g] = L[g,j] - a[j] - b[g] + m
# where L = log1p(UMI/sf) (sparse), a = colMeans(L), b = rowMeans(L), m = mean(L)
#
# irlba multiply functions (ctA is N x G, cells x genes):
#   ctA %*% v  = t(L) %*% v - (a - m) * sum(v) - dot(b, v) * 1_N
#   u %*% ctA  = (ctA^T %*% u)^T  [irlba then transposes this for t(ctA) %*% u]
#   ctA^T %*% u = L %*% u - dot(a, u) * 1_G - (b - m) * sum(u)
#
# Run from benchmark/ directory:
#   Rscript run_clr_gse184806.R

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(tibble)
})

set.seed(42)

RESULTS_DIR <- "output/benchmark_results"
DATA_DIR    <- "output/clr_local/data/consistency"

# ── Sparse CLR PCA ─────────────────────────────────────────────────────────────
clr_pca_sparse <- function(UMI, sf, n_pca) {
  stopifnot(inherits(UMI, "sparseMatrix"))
  G <- nrow(UMI)  # genes
  N <- ncol(UMI)  # cells

  # Step 1: L = log1p(UMI / sf) — sparse, same sparsity pattern as UMI
  L <- UMI
  col_idx <- rep(seq_len(N), diff(L@p))
  L@x <- log1p(L@x / sf[col_idx])

  # Step 2: Precompute row/col means for double-centering
  a <- Matrix::colMeans(L)  # per-cell mean, length N
  b <- Matrix::rowMeans(L)  # per-gene mean, length G
  m <- mean(a)               # grand mean (same as mean(b))

  a_c <- as.numeric(a - m)   # centered per-cell means
  b_c <- as.numeric(b - m)   # centered per-gene means
  a_n <- as.numeric(a)
  b_n <- as.numeric(b)

  ones_G <- rep(1.0, G)
  ones_N <- rep(1.0, N)

  # Step 3: Custom irlba multiply for ctA (N x G)
  # irlba calls: mult(A, v) for forward ctA %*% v
  #              mult(u, A) for reverse u %*% ctA (irlba transposes result)
  mult_fn <- function(x, y) {
    if (inherits(x, "clr_mat_")) {
      # Forward: ctA %*% y  (y is length-G vector or G x k matrix)
      if (is.matrix(y)) {
        # Apply column-by-column
        apply(y, 2L, function(v) {
          as.numeric(Matrix::t(L) %*% v - a_c * sum(v) - sum(b_n * v) * ones_N)
        })
      } else {
        v <- as.numeric(y)
        as.numeric(Matrix::t(L) %*% v - a_c * sum(v) - sum(b_n * v) * ones_N)
      }
    } else {
      # Reverse: x %*% ctA  (x is length-N vector or k x N matrix)
      # irlba transposes the result to get ctA^T %*% x^T
      if (is.matrix(x)) {
        apply(x, 1L, function(u) {
          as.numeric(L %*% u - sum(a_n * u) * ones_G - b_c * sum(u))
        })
      } else {
        u <- as.numeric(x)
        # Return as 1 x G row matrix (irlba transposes to get G x 1)
        matrix(as.numeric(L %*% u - sum(a_n * u) * ones_G - b_c * sum(u)), nrow = 1L)
      }
    }
  }

  # S3 class with dim methods so irlba gets the right dimensions
  A_fake <- structure(list(nrow_ = N, ncol_ = G), class = "clr_mat_")

  n <- min(n_pca, N - 1L, G - 1L)

  svd_res <- irlba(A = A_fake, nv = n, nu = n, mult = mult_fn)

  # PC scores (N x n_pca)
  svd_res$u %*% diag(svd_res$d, nrow = n, ncol = n)
}

# S3 dim methods for irlba's nrow(A) / ncol(A) calls
dim.clr_mat_  <- function(x) c(x$nrow_, x$ncol_)
nrow.clr_mat_ <- function(x) x$nrow_
ncol.clr_mat_ <- function(x) x$ncol_

make_knn <- function(emb, k) {
  k_use <- min(k, nrow(emb) - 1L)
  FNN::get.knnx(emb, emb, k = k_use + 1L)$nn.index[, -1L, drop = FALSE][, seq_len(k_use), drop = FALSE]
}

mean_overlap <- function(knn1, knn2) {
  mean(vapply(seq_len(nrow(knn1)), function(i) length(intersect(knn1[i,], knn2[i,])), 0.0))
}

size_factors <- function(UMI) { cs <- Matrix::colSums(UMI); cs / mean(cs) }

filter_mat <- function(UMI) UMI[Matrix::rowSums(UMI) > 0, Matrix::colSums(UMI) > 0, drop = FALSE]

KNN_VALS <- c(10, 50, 100)
PCA_VALS <- c(5, 10, 50)
SEEDS    <- 1:5

# ── Load GSE184806 ─────────────────────────────────────────────────────────────
message("Loading GSE184806 ...")
sce <- DropletUtils::read10xCounts(file.path(DATA_DIR, "GSE184806"))
UMI <- filter_mat(assay(sce, "counts"))
message("GSE184806: ", ncol(UMI), " cells x ", nrow(UMI), " genes")
gc()

# ── Quick sanity check on small subset ────────────────────────────────────────
message("Quick sanity check on 300-cell subset ...")
set.seed(999)
idx_t  <- sample(ncol(UMI), 300)
UMI_t  <- UMI[, idx_t]
sf_t   <- size_factors(UMI_t)
pca_t  <- clr_pca_sparse(UMI_t, sf_t, 5)
message("  Test PCA succeeded: ", nrow(pca_t), " cells x ", ncol(pca_t), " PCs")

# ── Run consistency ─────────────────────────────────────────────────────────
# Cache PCA per (seed, pca_dim) to avoid recomputing for each knn value
rows <- list()
for (seed in SEEDS) {
  set.seed(seed)
  idx1 <- sample(nrow(UMI), floor(nrow(UMI) / 2))
  idx2 <- setdiff(seq_len(nrow(UMI)), idx1)
  UMI1 <- UMI[idx1, ]; UMI2 <- UMI[idx2, ]
  sf1  <- size_factors(UMI1); sf2 <- size_factors(UMI2)

  for (pca_dim in PCA_VALS) {
    message(sprintf("  seed=%d pca=%d — computing PCA ...", seed, pca_dim))
    t_pca <- proc.time()
    pca1 <- clr_pca_sparse(UMI1, sf1, pca_dim)
    pca2 <- clr_pca_sparse(UMI2, sf2, pca_dim)
    elapsed_pca <- (proc.time() - t_pca)[["elapsed"]]
    message(sprintf("    PCA done in %.1fs", elapsed_pca))

    for (knn in KNN_VALS) {
      t0 <- proc.time()
      knn1 <- make_knn(pca1, knn)
      knn2 <- make_knn(pca2, knn)
      elapsed <- elapsed_pca + (proc.time() - t0)[["elapsed"]]
      ov <- mean_overlap(knn1, knn2)
      rows[[length(rows)+1]] <- tibble(
        mean_overlap      = ov,
        transformation_id = paste0("clr_GSE184806_s", seed, "_p", pca_dim, "_k", knn),
        dataset           = "GSE184806", seed = seed, pca_dim = pca_dim, knn = knn,
        transformation    = "clr", alpha = "FALSE",
        cputime_sec = elapsed, elapsed_sec = elapsed
      )
      message(sprintf("    knn=%d  overlap=%.3f", knn, ov))
    }
    gc()
  }
}

new_rows <- bind_rows(rows)
message("\nComputed ", nrow(new_rows), " CLR rows for GSE184806")

existing <- read_tsv(file.path(RESULTS_DIR, "consistency_results.tsv"),
                     show_col_types = FALSE) %>%
  mutate(alpha = as.character(alpha)) %>%
  filter(!(dataset == "GSE184806" & transformation == "clr"))

write_tsv(bind_rows(existing, new_rows),
          file.path(RESULTS_DIR, "consistency_results.tsv"))
message("Written to consistency_results.tsv")
message("=== DONE ===")
