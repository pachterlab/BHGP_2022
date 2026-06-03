#!/usr/bin/env Rscript
# Generate R-reference sctransform v2 (glmGamPoi) Pearson residuals for a raw
# count matrix, aligned to the FULL gene x cell frame (0 for genes/cells vst
# drops), written cells x genes to match the stored pysctransform output.
#
# Usage:
#   R_LD_LIBRARY_PATH=/usr/lib/R/lib:<buildtools>/lib \
#     Rscript gen_sctransform_v2.R <raw.mtx.gz> <out.csv.gz>
suppressMessages({ library(Matrix); library(sctransform); library(glmGamPoi) })

a <- commandArgs(trailingOnly = TRUE)
raw_fn <- a[[1]]; out_fn <- a[[2]]

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

# full is genes x cells (column-major). as.vector gives cell-blocks of genes,
# so np.fromfile(float32).reshape(n_cells, n_genes) yields a cells x genes
# matrix matching the stored sctransform.csv.gz orientation. Write float32 +
# a sidecar shape file (n_cells,n_genes). writeBin is ~seconds vs write.table.
con <- file(out_fn, "wb")
writeBin(as.vector(full), con, size = 4)
close(con)
writeLines(sprintf("%d,%d", ncol(M), nrow(M)), paste0(out_fn, ".shape"))
message("wrote ", out_fn, " (float32 binary, ", ncol(M), " cells x ", nrow(M), " genes)")
