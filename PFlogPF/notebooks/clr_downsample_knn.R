#!/usr/bin/env Rscript
# clr_downsample_knn.R
#
# Simple experiment: downsample ALL cells to target depth, build KNN,
# compare to KNN on full-depth data. No mixing, no integration.
#
# Usage: Rscript clr_downsample_knn.R <dataset> <target_depth>
#   dataset: "lung" (GSE150068) or "aorta" (GSE178765)

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(scuttle)
})

set.seed(42)

args         <- commandArgs(trailingOnly = TRUE)
DATASET      <- if (!is.na(args[1])) args[1] else "lung"
TARGET_DEPTH <- as.numeric(args[2])
if (is.na(TARGET_DEPTH)) TARGET_DEPTH <- 2000

N_HVG <- 2000

data_dirs <- list(
  lung  = "../benchmark/output/clr_local/data/consistency/GSM4522986",
  aorta = "../benchmark/output/clr_local/data/consistency/GSE178765"
)
DATA_DIR <- data_dirs[[DATASET]]

cat("Dataset:", DATASET, "\n")
cat("Target depth:", TARGET_DEPTH, "\n")

# ── Load ──────────────────────────────────────────────────────────────────────
mat   <- readMM(gzcon(file(file.path(DATA_DIR, "matrix.mtx.gz"), "rb")))
genes <- read.table(gzfile(file.path(DATA_DIR, "genes.tsv.gz")),
                    header = FALSE, stringsAsFactors = FALSE)$V1
cells <- read.table(gzfile(file.path(DATA_DIR, "barcodes.tsv.gz")),
                    header = FALSE, stringsAsFactors = FALSE)$V1
rownames(mat) <- genes; colnames(mat) <- cells
raw <- as(mat[Matrix::rowSums(mat) > 0, Matrix::colSums(mat) > 0], "dgCMatrix")
cat("Data:", nrow(raw), "genes x", ncol(raw), "cells\n")
cat("Median depth (full):", round(median(Matrix::colSums(raw))), "\n")

# ── Downsample ALL cells ──────────────────────────────────────────────────────
prop     <- TARGET_DEPTH / median(Matrix::colSums(raw))
raw_down <- scuttle::downsampleMatrix(raw, prop = prop, bycol = FALSE)
keep     <- Matrix::colSums(raw) > 0 & Matrix::colSums(raw_down) > 0
raw      <- raw[, keep]
raw_down <- raw_down[, keep]
cat("Median depth (downsampled):", round(median(Matrix::colSums(raw_down))), "\n")
cat("Cells retained:", ncol(raw), "\n")

# ── Transforms ────────────────────────────────────────────────────────────────
sf_vec <- function(UMI) { cs <- Matrix::colSums(UMI); cs / mean(cs) }

logp1_transform <- function(UMI) {
  as.matrix(log1p(Matrix::t(Matrix::t(UMI) / sf_vec(UMI))))
}

clr_transform <- function(UMI) {
  lpf <- logp1_transform(UMI)
  sweep(lpf, 2L, colMeans(lpf), "-")
}

select_hvg <- function(norm_mat, n = 2000) {
  gene_mean <- rowMeans(norm_mat)
  gene_var  <- apply(norm_mat, 1, var)
  lfit      <- loess(log(gene_var + 1e-10) ~ log(gene_mean + 1e-10), span = 0.3)
  order(log(gene_var + 1e-10) - fitted(lfit), decreasing = TRUE)[seq_len(n)]
}

run_pca <- function(mat, hvg_idx, n = 30) {
  prcomp_irlba(t(mat[hvg_idx, , drop = FALSE]), n = n,
               center = TRUE, scale. = FALSE)$x
}

make_knn <- function(emb, k = 50) {
  FNN::get.knnx(emb, emb, k = k + 1L)$nn.index[, -1L, drop = FALSE]
}

knn_overlap <- function(knn_a, knn_b, k = 50) {
  mean(vapply(seq_len(nrow(knn_a)), function(i)
    length(intersect(knn_a[i, ], knn_b[i, ])) / k, 0.0))
}

# ── Ground truth: full-depth logp1 KNN ───────────────────────────────────────
cat("\n── Ground truth (full depth, logp1, HVG) ──\n")
norm_full_logp1 <- logp1_transform(raw)
hvg_full        <- select_hvg(norm_full_logp1, n = N_HVG)
pca_full        <- run_pca(norm_full_logp1, hvg_full)
knn_full        <- make_knn(pca_full)
cat("  Built ground truth KNN on", ncol(raw), "cells\n")

# ── logp1 on downsampled ──────────────────────────────────────────────────────
cat("\n── logp1 (downsampled) ──\n")
norm_down_logp1 <- logp1_transform(raw_down)
hvg_down        <- select_hvg(norm_down_logp1, n = N_HVG)
pca_down_logp1  <- run_pca(norm_down_logp1, hvg_down)
knn_down_logp1  <- make_knn(pca_down_logp1)
ov_logp1        <- knn_overlap(knn_down_logp1, knn_full)
r_logp1         <- cor(pca_down_logp1[, 1], log10(Matrix::colSums(raw_down)))
cat("  KNN overlap:", round(ov_logp1, 4),
    "  |r(PC1,depth)|:", round(abs(r_logp1), 3), "\n")

# ── CLR on downsampled ────────────────────────────────────────────────────────
cat("\n── CLR (downsampled) ──\n")
norm_down_clr  <- clr_transform(raw_down)
pca_down_clr   <- run_pca(norm_down_clr, hvg_down)   # same HVG set as logp1
knn_down_clr   <- make_knn(pca_down_clr)
ov_clr         <- knn_overlap(knn_down_clr, knn_full)
r_clr          <- cor(pca_down_clr[, 1], log10(Matrix::colSums(raw_down)))
cat("  KNN overlap:", round(ov_clr, 4),
    "  |r(PC1,depth)|:", round(abs(r_clr), 3), "\n")

cat("\n=== Results (", DATASET, ", target depth:", TARGET_DEPTH, ") ===\n")
cat(sprintf("  %-10s  KNN=%.4f  |r|=%.3f\n", "logp1", ov_logp1, abs(r_logp1)))
cat(sprintf("  %-10s  KNN=%.4f  |r|=%.3f\n", "CLR",   ov_clr,   abs(r_clr)))
cat(sprintf("  CLR improvement: +%.1f%%\n", 100 * (ov_clr - ov_logp1) / ov_logp1))
