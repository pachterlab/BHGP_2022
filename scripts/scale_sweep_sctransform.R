#!/usr/bin/env Rscript
# Scale-invariance sweep for sctransform v2 (Satija recommendation).
#
# For a set of total-count scale factors s, scale the raw counts by s, round to
# integers, run sctransform::vst(vst.flavor="v2") (glmGamPoi backend), and report
# the RELATIVE Frobenius distance between the s-scaled residual matrix and the
# baseline (s = 1):  ||R(sM) - R(M)||_F / ||R(M)||_F.
#
# Residuals are aligned to the FULL gene set of the (sanitized) baseline; genes
# vst drops at a given scale contribute 0 residual, so every scale's matrix is
# the same shape and the Frobenius difference is well defined.
#
# Usage:
#   R_LD_LIBRARY_PATH=/usr/lib/R/lib:<buildtools>/lib \
#     Rscript scale_sweep_sctransform.R <raw.mtx.gz> <scales_csv> <out_tsv> [n_cells] [seed]
#
# <scales_csv> e.g. "0.1,0.2,0.5,1,2,5,10". <out_tsv> gets columns: scale, rel_frob.

suppressMessages({
  library(Matrix)
  library(sctransform)
  library(glmGamPoi)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: scale_sweep_sctransform.R raw.mtx.gz scales out.tsv [n_cells] [seed]")
raw_fn  <- args[[1]]
scales  <- as.numeric(strsplit(args[[2]], ",")[[1]])
out_tsv <- args[[3]]
n_cells <- if (length(args) >= 4) as.integer(args[[4]]) else 3000L
seed    <- if (length(args) >= 5) as.integer(args[[5]]) else 0L

# raw.mtx.gz is cells x genes (MatrixMarket). vst wants genes x cells.
message("reading ", raw_fn)
M <- as(readMM(raw_fn), "CsparseMatrix")          # cells x genes
M <- t(M)                                          # genes x cells

# Subsample cells (columns) for tractability; fixed across all scales.
set.seed(seed)
if (ncol(M) > n_cells) {
  idx <- sort(sample(ncol(M), n_cells))
  M <- M[, idx]
}
# Sanitize once: drop all-zero genes/cells so the baseline gene set is fixed.
M <- M[Matrix::rowSums(M) > 0, , drop = FALSE]
M <- M[, Matrix::colSums(M) > 0, drop = FALSE]
rownames(M) <- paste0("g", seq_len(nrow(M)))
colnames(M) <- paste0("c", seq_len(ncol(M)))
all_genes <- rownames(M)
all_cells <- colnames(M)
message(sprintf("matrix: %d genes x %d cells", nrow(M), ncol(M)))

# Residuals for one scale, reindexed to the fixed gene x cell frame. Cells/genes
# that drop to all-zero after scaling+rounding (or that vst filters) get residual
# 0 -- a zeroed cell carries no signal, so 0 keeps every scale the same shape.
resid_full <- function(scaled) {
  keep_g <- Matrix::rowSums(scaled) > 0
  keep_c <- Matrix::colSums(scaled) > 0
  sub <- scaled[keep_g, keep_c, drop = FALSE]
  v <- vst(sub, vst.flavor = "v2", verbosity = 0,
           method = "glmGamPoi_offset", return_cell_attr = FALSE)
  y <- v$y                                          # (kept genes) x (kept cells)
  full <- matrix(0.0, nrow = length(all_genes), ncol = length(all_cells),
                 dimnames = list(all_genes, all_cells))
  full[rownames(y), colnames(y)] <- y
  full
}

scale_counts <- function(s) {
  sc <- M * s
  sc@x <- round(sc@x)
  sc <- drop0(sc)
  sc
}

message("baseline (s=1)")
R1 <- resid_full(scale_counts(1))
base_norm <- sqrt(sum(R1^2))

res <- data.frame(scale = numeric(0), rel_frob = numeric(0))
for (s in scales) {
  message(sprintf("scale %.4g", s))
  Rs <- if (s == 1) R1 else resid_full(scale_counts(s))
  d  <- sqrt(sum((Rs - R1)^2)) / base_norm
  res <- rbind(res, data.frame(scale = s, rel_frob = d))
  cat(sprintf("  s=%.4g  rel_frob=%.5f\n", s, d))
}

write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", out_tsv)
