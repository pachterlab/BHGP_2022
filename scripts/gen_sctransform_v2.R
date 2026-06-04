#!/usr/bin/env Rscript
# Generate R-reference sctransform v2 (glmGamPoi) Pearson residuals for a raw
# count matrix, aligned to the FULL gene x cell frame (0 for genes/cells vst
# drops), written cells x genes to match the stored pysctransform output.
#
# Output format is by extension: .csv.gz -> data.table::fwrite (fast, gzip,
# cells x genes, no header/rownames; matches stored sctransform.csv.gz);
# .bin -> float32 binary + .shape sidecar.
#
# Usage:
#   R_LD_LIBRARY_PATH=/usr/lib/R/lib:<buildtools>/lib \
#     Rscript gen_sctransform_v2.R <raw.mtx.gz> <out.csv.gz|out.bin> [n_threads]
suppressMessages({ library(Matrix); library(sctransform); library(glmGamPoi) })

# sctransform parallelizes parameter estimation via 'future'; we already
# parallelize across datasets, so force sequential (avoids oversubscription and
# the default 500 MiB future.globals.maxSize error on larger matrices).
suppressMessages(if (requireNamespace("future", quietly = TRUE)) {
  options(future.globals.maxSize = 16 * 1024^3)
  future::plan("sequential")
})

a <- commandArgs(trailingOnly = TRUE)
raw_fn <- a[[1]]; out_fn <- a[[2]]
nthreads <- if (length(a) >= 3) as.integer(a[[3]]) else 2L

M <- t(as(readMM(raw_fn), "CsparseMatrix"))          # genes x cells
rownames(M) <- paste0("g", seq_len(nrow(M)))
colnames(M) <- paste0("c", seq_len(ncol(M)))
message(sprintf("raw: %d genes x %d cells", nrow(M), ncol(M)))

keep_g <- Matrix::rowSums(M) > 0
keep_c <- Matrix::colSums(M) > 0
v <- vst(M[keep_g, keep_c, drop = FALSE], vst.flavor = "v2",
         method = "glmGamPoi_offset", verbosity = 0)
y <- v$y                                             # kept genes x kept cells

full <- matrix(0.0, nrow = nrow(M), ncol = ncol(M),
               dimnames = list(rownames(M), colnames(M)))
full[rownames(y), colnames(y)] <- y
message(sprintf("residuals placed; kept %d genes, %d cells", nrow(y), ncol(y)))

if (grepl("\\.csv\\.gz$", out_fn)) {
  # cells x genes, no header / no row names; matches stored sctransform.csv.gz.
  suppressMessages(library(data.table))
  data.table::setDTthreads(nthreads)
  data.table::fwrite(data.table::as.data.table(t(full)), out_fn,
                     col.names = FALSE, compress = "gzip")
  message("wrote ", out_fn, " (csv.gz, ", ncol(M), " cells x ", nrow(M), " genes)")
} else {
  # full is genes x cells (column-major). as.vector gives cell-blocks of genes,
  # so np.fromfile(float32).reshape(n_cells, n_genes) yields a cells x genes
  # matrix matching the sctransform.csv.gz orientation.
  con <- file(out_fn, "wb")
  writeBin(as.vector(full), con, size = 4)
  close(con)
  writeLines(sprintf("%d,%d", ncol(M), nrow(M)), paste0(out_fn, ".shape"))
  message("wrote ", out_fn, " (float32 binary, ", ncol(M), " cells x ", nrow(M), " genes)")
}
