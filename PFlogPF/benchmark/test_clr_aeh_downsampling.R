#!/usr/bin/env Rscript
# Sanity test: re-run CLR downsampling on the mcSCRB dataset using AE&H's
# downsampling methodology and compare to the existing per-cell rmultinom
# result currently in downsampling_results.tsv.
#
# AE&H method:  scuttle::downsampleMatrix(UMI, prop = 5000/median(colsums), bycol = FALSE)
#               Each matrix entry is independently downsampled: rbinom(x_ij, prop).
# Old CLR method: per-cell rmultinom(size = sum(cell) * 0.1, prob = cell + 1e-8)
#
# Single dataset (mcSCRB), 5 seeds, pca_dim = 10, knn = 50 — same parameter
# slice used in the main panel of the benchmark figure.

suppressPackageStartupMessages({
  library(Matrix); library(irlba); library(FNN)
})

set.seed(42)
DATA <- "output/clr_local/data/downsampling/GSE103568_JM8_UMIcounts.txt.gz"
URL  <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103568&format=file&file=GSE103568%5FJM8%5FUMIcounts%2Etxt%2Egz"

dir.create(dirname(DATA), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(DATA)) {
  message("Downloading mcSCRB ...")
  download.file(URL, DATA, mode = "wb", quiet = TRUE)
}

UMI_full <- as.matrix(read.delim(DATA))
expressed_genes <- rowSums(UMI_full) > 0
expressed_cells <- colSums(UMI_full) > 0
UMI_full <- UMI_full[expressed_genes, expressed_cells]
storage.mode(UMI_full) <- "integer"
message(sprintf("mcSCRB: %d cells x %d genes; median(colsums) = %.0f, mean = %.0f",
                ncol(UMI_full), nrow(UMI_full),
                median(colSums(UMI_full)), mean(colSums(UMI_full))))

# scuttle::downsampleMatrix(..., bycol=FALSE) is equivalent to: for every entry
# x_ij, draw rbinom(1, x_ij, prop). Implemented manually here so we don't depend
# on Bioconductor's scuttle install (renv mismatch makes that painful).
downsample_matrix <- function(M, prop) {
  if (!is.matrix(M)) M <- as.matrix(M)
  out <- matrix(rbinom(length(M), as.vector(M), prop), nrow = nrow(M))
  dimnames(out) <- dimnames(M)
  storage.mode(out) <- "integer"
  out
}

# Per-cell rmultinom downsampling (the original CLR pipeline's approach).
downsample_per_cell <- function(M, fraction = 0.1) {
  out <- apply(M, 2, function(cell) {
    n <- max(1L, round(sum(cell) * fraction))
    as.integer(rmultinom(1L, size = n, prob = cell + 1e-8))
  })
  rownames(out) <- rownames(M)
  out
}

# CLR transform + PCA + kNN, returning the kNN index matrix.
clr_knn <- function(UMI, pca_dim = 10, knn = 50) {
  sf <- colSums(UMI); sf <- sf / mean(sf)
  logpf <- log1p(sweep(UMI, 2L, sf, "/"))
  X <- t(sweep(logpf, 2L, colMeans(logpf), "-"))
  k <- min(pca_dim, nrow(X) - 1L, ncol(X) - 1L)
  pc <- prcomp_irlba(X, n = k, center = TRUE, scale. = FALSE)$x
  FNN::get.knnx(pc, pc, k = knn + 1L)$nn.index[, -1L, drop = FALSE]
}

mean_knn_overlap <- function(a, b) {
  mean(vapply(seq_len(nrow(a)),
              function(i) length(intersect(a[i,], b[i,])), 0.0))
}

# AE&H prop for this dataset.
prop_aeh <- 5000 / median(colSums(UMI_full))
message(sprintf("AE&H prop = 5000/median = %.4f  (vs per-cell fraction = 0.1)", prop_aeh))

results <- data.frame(seed = integer(), method = character(), overlap = numeric())
for (seed in 1:5) {
  set.seed(seed)
  # AE&H method
  UMI_aeh <- downsample_matrix(UMI_full, prop_aeh)
  UMI_aeh <- UMI_aeh[rowSums(UMI_aeh) > 0, colSums(UMI_aeh) > 0]
  cg <- intersect(rownames(UMI_full), rownames(UMI_aeh))
  k_full_a <- clr_knn(UMI_full[cg, ], pca_dim = 10, knn = 50)
  k_red_a  <- clr_knn(UMI_aeh[cg, ],  pca_dim = 10, knn = 50)
  ov_aeh <- mean_knn_overlap(k_full_a, k_red_a)

  # Per-cell rmultinom method
  set.seed(seed)
  UMI_pc <- downsample_per_cell(UMI_full, fraction = 0.1)
  UMI_pc <- UMI_pc[rowSums(UMI_pc) > 0, colSums(UMI_pc) > 0]
  cg <- intersect(rownames(UMI_full), rownames(UMI_pc))
  k_full_p <- clr_knn(UMI_full[cg, ], pca_dim = 10, knn = 50)
  k_red_p  <- clr_knn(UMI_pc[cg, ],   pca_dim = 10, knn = 50)
  ov_pc <- mean_knn_overlap(k_full_p, k_red_p)

  message(sprintf("seed=%d: AE&H overlap=%.3f   per-cell overlap=%.3f", seed, ov_aeh, ov_pc))
  results <- rbind(results,
    data.frame(seed = seed, method = "aeh_downsampleMatrix",   overlap = ov_aeh),
    data.frame(seed = seed, method = "percell_rmultinom_0.1", overlap = ov_pc))
}

cat("\n=== Summary (knn=50, pca_dim=10, mcSCRB) ===\n")
agg <- aggregate(overlap ~ method, data = results, function(x)
  c(min = min(x), mean = mean(x), max = max(x)))
print(agg)
cat("\nExpected: other transformations land near overlap ~6-7 (12-14% of knn).\n")
cat("If AE&H-method CLR is also in that range, the CLR advantage was a methodology artifact.\n")
cat("If AE&H-method CLR stays in the 30s, the CLR advantage is real.\n")
