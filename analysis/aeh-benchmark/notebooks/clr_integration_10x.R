#!/usr/bin/env Rscript
# clr_integration_10x.R
#
# Same experiment as clr_harmony_integration.R / clr_seurat_integration.R
# but on a large 10x dataset (GSE178765, Mouse Aorta, ~10k cells).
# Uses RPCA for Seurat (faster and more scalable than CCA at this scale).
#
# Usage: Rscript clr_integration_10x.R <target_depth>

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(dplyr)
  library(harmony)
  library(scuttle)
})

set.seed(42)
options(future.globals.maxSize = Inf)  # no limit for Seurat parallel

DATA_DIR <- "../benchmark/output/clr_local/data/consistency/GSM4522986"
OUT_DIR  <- "../output"
dir.create(OUT_DIR, showWarnings = FALSE)

# ── 1. Load 10x MTX data ──────────────────────────────────────────────────────
cat("Loading GSM4522986 (GSE150068 Lung Epithelium 10x)...\n")
mat   <- readMM(gzcon(file(file.path(DATA_DIR, "matrix.mtx.gz"), "rb")))
genes <- read.table(gzfile(file.path(DATA_DIR, "genes.tsv.gz")),
                    header = FALSE, stringsAsFactors = FALSE)$V1
cells <- read.table(gzfile(file.path(DATA_DIR, "barcodes.tsv.gz")),
                    header = FALSE, stringsAsFactors = FALSE)$V1
rownames(mat) <- genes
colnames(mat) <- cells
raw <- mat[Matrix::rowSums(mat) > 0, Matrix::colSums(mat) > 0]
raw <- as(raw, "dgCMatrix")
cat("  Data:", nrow(raw), "genes x", ncol(raw), "cells\n")
cat("  Median depth:", round(median(Matrix::colSums(raw))), "reads/cell\n")

# ── 2. Downsample ─────────────────────────────────────────────────────────────
TARGET_DEPTH <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
if (is.na(TARGET_DEPTH)) TARGET_DEPTH <- 2000
cat("Target depth:", TARGET_DEPTH, "reads/cell\n")
prop <- TARGET_DEPTH / median(Matrix::colSums(raw))
cat("  Proportion:", round(prop, 3), "\n")
raw_down <- scuttle::downsampleMatrix(raw, prop = prop, bycol = FALSE)
keep <- Matrix::colSums(raw) > 0 & Matrix::colSums(raw_down) > 0
raw      <- raw[, keep]
raw_down <- raw_down[, keep]
cat("  Cells retained:", ncol(raw), "\n")
cat("  Median depth — full:", round(median(Matrix::colSums(raw))),
    "  downsampled:", round(median(Matrix::colSums(raw_down))), "\n")

# ── 3. Split into two batches ─────────────────────────────────────────────────
n_cells  <- ncol(raw)
half     <- floor(n_cells / 2)
idx_full <- sample(n_cells, half)
idx_down <- setdiff(seq_len(n_cells), idx_full)
orig_idx <- c(idx_full, idx_down)

UMI_gt  <- raw
UMI_mix <- cbind(raw[, idx_full], raw_down[, idx_down])
batch_labels <- c(rep("full", length(idx_full)), rep("downsampled", length(idx_down)))
colnames(UMI_mix) <- paste0("cell", seq_len(ncol(UMI_mix)))

cat("Mixed dataset:", ncol(UMI_mix), "cells (",
    sum(batch_labels == "full"), "full +",
    sum(batch_labels == "downsampled"), "downsampled)\n")

# ── 4. Helper functions ───────────────────────────────────────────────────────
sf_vec <- function(UMI) { cs <- Matrix::colSums(UMI); cs / mean(cs) }

logp1_transform <- function(UMI) {
  sf  <- sf_vec(UMI)
  log1p(Matrix::t(Matrix::t(UMI) / sf))
}

clr_transform <- function(UMI) {
  sf  <- sf_vec(UMI)
  lpf <- log1p(Matrix::t(Matrix::t(UMI) / sf))
  lpf - Matrix::rowMeans(lpf)   # subtract per-cell mean (lpf is genes x cells; rowMeans = per-gene, so need colMeans)
}

# Fix: CLR subtracts per-CELL mean, i.e. colMeans of the genes×cells matrix
clr_transform <- function(UMI) {
  sf  <- sf_vec(UMI)
  lpf <- as.matrix(log1p(Matrix::t(Matrix::t(UMI) / sf)))
  sweep(lpf, 2L, colMeans(lpf), "-")
}

logp1_transform <- function(UMI) {
  sf  <- sf_vec(UMI)
  as.matrix(log1p(Matrix::t(Matrix::t(UMI) / sf)))
}

select_hvg <- function(norm_mat, n = 2000) {
  # Standard HVG selection: mean and variance of each gene across cells,
  # fit loess mean-variance trend, pick top n by residual variance.
  gene_mean <- rowMeans(norm_mat)
  gene_var  <- apply(norm_mat, 1, var)
  # log-transform for stable loess fit (add small offset to avoid log(0))
  lfit <- loess(log(gene_var + 1e-10) ~ log(gene_mean + 1e-10), span = 0.3)
  resid_var <- log(gene_var + 1e-10) - fitted(lfit)
  order(resid_var, decreasing = TRUE)[seq_len(n)]
}

run_pca_mat <- function(mat_genes_x_cells, n = 30, hvg_idx = NULL) {
  if (!is.null(hvg_idx)) mat_genes_x_cells <- mat_genes_x_cells[hvg_idx, , drop = FALSE]
  prcomp_irlba(t(mat_genes_x_cells), n = n, center = TRUE, scale. = FALSE)$x
}

make_knn <- function(embedding, k = 50) {
  nn <- FNN::get.knnx(embedding, embedding, k = k + 1L)$nn.index
  nn[, -1L, drop = FALSE]
}

overlap_with_gt <- function(knn_mixed, knn_gt, orig_idx, k = 50) {
  n <- nrow(knn_mixed)
  mean(vapply(seq_len(n), function(i) {
    gt_n  <- knn_gt[orig_idx[i], ]
    mix_n <- orig_idx[knn_mixed[i, ]]
    length(intersect(gt_n, mix_n)) / k
  }, 0.0))
}

run_harmony_int <- function(pca, batch) {
  harmony::HarmonyMatrix(pca, meta_data = data.frame(batch = batch),
                         vars_use = "batch", do_pca = FALSE, verbose = FALSE)
}

run_seurat_rpca <- function(norm_mat, raw_counts_mat, batch,
                             n_pcs = 30, n_features = 2000) {
  batches  <- unique(batch)

  # First pass: find shared HVGs across batches to reduce matrix size
  hvg_list <- lapply(batches, function(b) {
    idx  <- which(batch == b)
    seu_tmp <- CreateSeuratObject(counts = raw_counts_mat[, idx, drop = FALSE])
    seu_tmp[["RNA"]]@data <- as(norm_mat[, idx, drop = FALSE], "dgCMatrix")
    seu_tmp <- FindVariableFeatures(seu_tmp, selection.method = "vst",
                                    nfeatures = n_features, verbose = FALSE)
    VariableFeatures(seu_tmp)
  })
  shared_hvg <- Reduce(intersect, hvg_list)
  if (length(shared_hvg) < n_pcs * 2) shared_hvg <- unique(unlist(hvg_list))[seq_len(n_features)]
  shared_hvg <- shared_hvg[seq_len(min(n_features, length(shared_hvg)))]

  # Second pass: build Seurat objects restricted to HVGs
  seu_list <- lapply(batches, function(b) {
    idx      <- which(batch == b)
    counts_b <- raw_counts_mat[shared_hvg, idx, drop = FALSE]
    norm_b   <- as(norm_mat[shared_hvg, idx, drop = FALSE], "dgCMatrix")
    seu      <- CreateSeuratObject(counts = counts_b)
    seu[["RNA"]]@data <- norm_b
    VariableFeatures(seu) <- shared_hvg
    seu <- ScaleData(seu, features = shared_hvg, verbose = FALSE)
    seu <- RunPCA(seu, features = shared_hvg, npcs = n_pcs, verbose = FALSE)
    seu
  })
  names(seu_list) <- batches

  anchors    <- FindIntegrationAnchors(
    object.list     = seu_list,
    dims            = seq_len(n_pcs),
    anchor.features = shared_hvg,
    reduction       = "rpca",
    verbose         = FALSE
  )
  integrated <- IntegrateData(anchorset = anchors, dims = seq_len(n_pcs),
                               verbose = FALSE)
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = n_pcs, verbose = FALSE)

  emb        <- Embeddings(integrated, "pca")
  cell_order <- c(colnames(seu_list[["full"]]), colnames(seu_list[["downsampled"]]))
  emb[cell_order, , drop = FALSE]
}

# ── 5. Ground truth KNN ───────────────────────────────────────────────────────
cat("\n── Ground truth KNN (all", ncol(UMI_gt), "cells, full depth) ──\n")
norm_gt    <- logp1_transform(UMI_gt)
hvg_gt     <- select_hvg(norm_gt, n = 2000)
cat("  HVGs selected:", length(hvg_gt), "\n")
pca_gt <- run_pca_mat(norm_gt, n = 30, hvg_idx = hvg_gt)
knn_gt <- make_knn(pca_gt, k = 50)

# Select HVGs from the mixed dataset logp1 (shared gene set for both methods)
norm_logp1 <- logp1_transform(UMI_mix)
hvg_mix    <- select_hvg(norm_logp1, n = 2000)
cat("  HVGs selected (mixed):", length(hvg_mix), "\n")

# ── 6. logp1 — no integration ─────────────────────────────────────────────────
cat("\n── logp1 (no integration) ──\n")
pca_logp1         <- run_pca_mat(norm_logp1, n = 30, hvg_idx = hvg_mix)
knn_logp1_nointe  <- make_knn(pca_logp1, k = 50)
ov_logp1_nointe   <- overlap_with_gt(knn_logp1_nointe, knn_gt, orig_idx)
r_logp1_nointe    <- cor(pca_logp1[, 1], log10(Matrix::colSums(UMI_mix)))
cat("  KNN overlap:", round(ov_logp1_nointe, 4),
    "  |r(PC1,depth)|:", round(abs(r_logp1_nointe), 3), "\n")

# ── 7. logp1 + Harmony ────────────────────────────────────────────────────────
cat("\n── logp1 + Harmony ──\n")
harm_logp1       <- run_harmony_int(pca_logp1, batch_labels)
knn_logp1_h      <- make_knn(harm_logp1, k = 50)
ov_logp1_h       <- overlap_with_gt(knn_logp1_h, knn_gt, orig_idx)
r_logp1_h        <- cor(harm_logp1[, 1], log10(Matrix::colSums(UMI_mix)))
cat("  KNN overlap:", round(ov_logp1_h, 4),
    "  |r(PC1,depth)|:", round(abs(r_logp1_h), 3), "\n")

# ── 9. CLR — no integration ───────────────────────────────────────────────────
cat("\n── CLR (no integration) ──\n")
norm_clr         <- clr_transform(UMI_mix)
pca_clr          <- run_pca_mat(norm_clr, n = 30, hvg_idx = hvg_mix)
knn_clr_nointe   <- make_knn(pca_clr, k = 50)
ov_clr_nointe    <- overlap_with_gt(knn_clr_nointe, knn_gt, orig_idx)
r_clr_nointe     <- cor(pca_clr[, 1], log10(Matrix::colSums(UMI_mix)))
cat("  KNN overlap:", round(ov_clr_nointe, 4),
    "  |r(PC1,depth)|:", round(abs(r_clr_nointe), 3), "\n")

# ── 10. CLR + Harmony ─────────────────────────────────────────────────────────
cat("\n── CLR + Harmony ──\n")
harm_clr       <- run_harmony_int(pca_clr, batch_labels)
knn_clr_h      <- make_knn(harm_clr, k = 50)
ov_clr_h       <- overlap_with_gt(knn_clr_h, knn_gt, orig_idx)
r_clr_h        <- cor(harm_clr[, 1], log10(Matrix::colSums(UMI_mix)))
cat("  KNN overlap:", round(ov_clr_h, 4),
    "  |r(PC1,depth)|:", round(abs(r_clr_h), 3), "\n")

# ── 11. Summary ───────────────────────────────────────────────────────────────
cat("\n=== Results (target depth:", TARGET_DEPTH, ") ===\n")
cat(sprintf("  %-30s  KNN=%.4f  |r|=%.3f\n", "logp1 (no integration)", ov_logp1_nointe, abs(r_logp1_nointe)))
cat(sprintf("  %-30s  KNN=%.4f  |r|=%.3f\n", "logp1 + Harmony",        ov_logp1_h,      abs(r_logp1_h)))
cat(sprintf("  %-30s  KNN=%.4f  |r|=%.3f\n", "CLR (no integration)",   ov_clr_nointe,   abs(r_clr_nointe)))
cat(sprintf("  %-30s  KNN=%.4f  |r|=%.3f\n", "CLR + Harmony",          ov_clr_h,        abs(r_clr_h)))
