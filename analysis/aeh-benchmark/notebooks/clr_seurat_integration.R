#!/usr/bin/env Rscript
# clr_seurat_integration.R
#
# Same experiment as clr_harmony_integration.R but using Seurat CCA integration
# instead of Harmony.
#
# Setup:
#   - Take mcSCRB dataset (full depth)
#   - Downsample to TARGET_DEPTH reads/cell
#   - Mix: 50% full-depth + 50% downsampled (same cells)
#   - Integrate with Seurat (FindIntegrationAnchors + IntegrateData)
#   - Measure KNN overlap vs ground truth = KNN on ALL full-depth cells

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(Seurat)
  library(scuttle)
})

set.seed(42)

DATA_FILE <- "../benchmark/output/clr_local/data/downsampling/GSE103568_JM8_UMIcounts.txt.gz"
OUT_DIR   <- "../output"
dir.create(OUT_DIR, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading mcSCRB data...\n")
raw <- as.matrix(read.delim(DATA_FILE))
raw <- raw[rowSums(raw) > 0, colSums(raw) > 0]
cat("  Full data:", nrow(raw), "genes x", ncol(raw), "cells\n")
cat("  Median depth:", median(colSums(raw)), "reads/cell\n")

# ── 2. Downsample ─────────────────────────────────────────────────────────────
TARGET_DEPTH <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
if (is.na(TARGET_DEPTH)) TARGET_DEPTH <- 2000
cat("Downsampling to target depth:", TARGET_DEPTH, "reads/cell\n")
prop <- TARGET_DEPTH / median(colSums(raw))
cat("  Proportion:", round(prop, 3), "\n")
raw_down <- as.matrix(scuttle::downsampleMatrix(raw, prop = prop, bycol = FALSE))
keep <- colSums(raw) > 0 & colSums(raw_down) > 0
raw      <- raw[, keep]
raw_down <- raw_down[, keep]
cat("  Cells retained:", ncol(raw), "\n")
cat("  Median depth — full:", round(median(colSums(raw))),
    "  downsampled:", round(median(colSums(raw_down))), "\n")

# ── 3. Split into two batches ─────────────────────────────────────────────────
n_cells  <- ncol(raw)
half     <- floor(n_cells / 2)
idx_full <- sample(n_cells, half)
idx_down <- setdiff(seq_len(n_cells), idx_full)
orig_idx <- c(idx_full, idx_down)

UMI_gt  <- raw
UMI_mix <- cbind(raw[, idx_full], raw_down[, idx_down])
batch_labels <- c(rep("full", length(idx_full)), rep("downsampled", length(idx_down)))

cat("Mixed dataset:", ncol(UMI_mix), "cells (",
    sum(batch_labels=="full"), "full +",
    sum(batch_labels=="downsampled"), "downsampled)\n")

# ── 4. Helper functions ───────────────────────────────────────────────────────
sf_vec <- function(UMI) { cs <- colSums(UMI); cs / mean(cs) }

logp1_transform <- function(UMI) {
  sf <- sf_vec(UMI)
  log1p(sweep(UMI, 2L, sf, "/"))
}

clr_transform <- function(UMI) {
  sf  <- sf_vec(UMI)
  lpf <- log1p(sweep(UMI, 2L, sf, "/"))
  sweep(lpf, 2L, colMeans(lpf), "-")
}

run_pca_mat <- function(mat_genes_x_cells, n = 30) {
  prcomp_irlba(t(mat_genes_x_cells), n = n, center = TRUE, scale. = FALSE)$x
}

make_knn <- function(embedding, k = 50) {
  nn <- FNN::get.knnx(embedding, embedding, k = k + 1L)$nn.index
  nn[, -1L, drop = FALSE]
}

overlap_with_gt <- function(knn_mixed, knn_gt, orig_idx, k = 50) {
  n <- nrow(knn_mixed)
  mean(vapply(seq_len(n), function(i) {
    gt_neighbors  <- knn_gt[orig_idx[i], ]
    mix_neighbors <- orig_idx[knn_mixed[i, ]]
    length(intersect(gt_neighbors, mix_neighbors)) / k
  }, 0.0))
}

# ── 5. Seurat integration helper ──────────────────────────────────────────────
# Takes a pre-computed normalized matrix (genes x cells) and the batch labels,
# runs FindIntegrationAnchors + IntegrateData, returns corrected PCA embedding.
run_seurat_integration <- function(norm_mat, raw_counts_mat, batch,
                                    n_pcs = 30, n_features = 2000) {
  # We supply custom normalized data (logp1 or CLR) in the data slot.
  # Raw UMI counts go in the counts slot so FindVariableFeatures VST works.
  batches <- unique(batch)

  seu_list <- lapply(batches, function(b) {
    idx <- which(batch == b)
    counts_b <- Matrix::Matrix(raw_counts_mat[, idx, drop = FALSE], sparse = TRUE)
    norm_b   <- as(norm_mat[, idx, drop = FALSE], "dgCMatrix")
    seu <- CreateSeuratObject(counts = counts_b)
    seu[["RNA"]]@data <- norm_b
    seu <- FindVariableFeatures(seu, selection.method = "vst",
                                nfeatures = n_features, verbose = FALSE)
    seu
  })
  names(seu_list) <- batches

  anchors <- FindIntegrationAnchors(
    object.list      = seu_list,
    dims             = seq_len(n_pcs),
    anchor.features  = n_features,
    reduction        = "cca",
    verbose          = FALSE
  )
  integrated <- IntegrateData(anchorset = anchors, dims = seq_len(n_pcs),
                               verbose = FALSE)
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = n_pcs, verbose = FALSE)

  # Return PCA embedding in original cell order (Seurat may reorder)
  emb <- Embeddings(integrated, "pca")
  # Reconstruct original order: full-batch cells first, then downsampled
  cell_order <- c(colnames(seu_list[["full"]]), colnames(seu_list[["downsampled"]]))
  emb[cell_order, , drop = FALSE]
}

# ── 6. Ground truth KNN ───────────────────────────────────────────────────────
cat("\n── Ground truth KNN (all", ncol(UMI_gt), "cells, full depth) ──\n")
pca_gt <- run_pca_mat(logp1_transform(UMI_gt), n = 30)
knn_gt <- make_knn(pca_gt, k = 50)

# ── 7. logp1 — no integration ─────────────────────────────────────────────────
cat("\n── logp1 (no integration) ──\n")
pca_logp1        <- run_pca_mat(logp1_transform(UMI_mix), n = 30)
knn_logp1_nointe <- make_knn(pca_logp1, k = 50)
ov_logp1_nointe  <- overlap_with_gt(knn_logp1_nointe, knn_gt, orig_idx)
r_logp1_nointe   <- cor(pca_logp1[, 1], log10(colSums(UMI_mix)))
cat("  KNN overlap:", round(ov_logp1_nointe, 4),
    "  |r(PC1,depth)|:", round(abs(r_logp1_nointe), 3), "\n")

# ── 8. logp1 + Seurat integration ─────────────────────────────────────────────
cat("\n── logp1 + Seurat ──\n")
norm_logp1   <- logp1_transform(UMI_mix)
colnames(norm_logp1) <- paste0("cell", seq_len(ncol(norm_logp1)))
raw_mix_named <- UMI_mix; colnames(raw_mix_named) <- paste0("cell", seq_len(ncol(UMI_mix)))
seu_logp1    <- run_seurat_integration(norm_logp1, raw_mix_named, batch_labels)
knn_logp1_s  <- make_knn(seu_logp1, k = 50)
ov_logp1_s   <- overlap_with_gt(knn_logp1_s, knn_gt, orig_idx)
r_logp1_s    <- cor(seu_logp1[, 1], log10(colSums(UMI_mix)))
cat("  KNN overlap:", round(ov_logp1_s, 4),
    "  |r(PC1,depth)|:", round(abs(r_logp1_s), 3), "\n")

# ── 9. CLR — no integration ───────────────────────────────────────────────────
cat("\n── CLR (no integration) ──\n")
pca_clr        <- run_pca_mat(clr_transform(UMI_mix), n = 30)
knn_clr_nointe <- make_knn(pca_clr, k = 50)
ov_clr_nointe  <- overlap_with_gt(knn_clr_nointe, knn_gt, orig_idx)
r_clr_nointe   <- cor(pca_clr[, 1], log10(colSums(UMI_mix)))
cat("  KNN overlap:", round(ov_clr_nointe, 4),
    "  |r(PC1,depth)|:", round(abs(r_clr_nointe), 3), "\n")

# ── 10. CLR + Seurat integration ──────────────────────────────────────────────
cat("\n── CLR + Seurat ──\n")
norm_clr   <- clr_transform(UMI_mix)
colnames(norm_clr) <- paste0("cell", seq_len(ncol(norm_clr)))
seu_clr    <- run_seurat_integration(norm_clr, raw_mix_named, batch_labels)
knn_clr_s  <- make_knn(seu_clr, k = 50)
ov_clr_s   <- overlap_with_gt(knn_clr_s, knn_gt, orig_idx)
r_clr_s    <- cor(seu_clr[, 1], log10(colSums(UMI_mix)))
cat("  KNN overlap:", round(ov_clr_s, 4),
    "  |r(PC1,depth)|:", round(abs(r_clr_s), 3), "\n")

# ── 11. Summary ───────────────────────────────────────────────────────────────
cat("\n=== Results (target depth:", TARGET_DEPTH, ") ===\n")
cat(sprintf("  %-28s  KNN=%.4f  |r|=%.3f\n",
            "logp1 (no Seurat)", ov_logp1_nointe, abs(r_logp1_nointe)))
cat(sprintf("  %-28s  KNN=%.4f  |r|=%.3f\n",
            "logp1 + Seurat", ov_logp1_s, abs(r_logp1_s)))
cat(sprintf("  %-28s  KNN=%.4f  |r|=%.3f\n",
            "CLR (no Seurat)", ov_clr_nointe, abs(r_clr_nointe)))
cat(sprintf("  %-28s  KNN=%.4f  |r|=%.3f\n",
            "CLR + Seurat", ov_clr_s, abs(r_clr_s)))
