#!/usr/bin/env Rscript
# Pilot: empirical CV(q) for sctransform under two gene-set choices, using the
# pipeline's exact vst call (vst.flavor="v2", method="glmGamPoi_offset").
#   full-fit          : vst on full subset matrix; CV over ALL kept genes (== stored cov_gene)
#   full-fit-filtered : same residuals, CV over filtered genes (mu>0.05 & det>=1%)
#   filtered-fit      : vst RE-RUN on the filtered-gene matrix; CV over its genes
# args: <dataset> <do_full 0|1>
suppressMessages({ library(Matrix); library(sctransform); library(glmGamPoi) })
suppressMessages(if (requireNamespace("future", quietly=TRUE)) {
  options(future.globals.maxSize=16*1024^3); future::plan("sequential") })
a <- commandArgs(trailingOnly=TRUE); ds <- a[[1]]; do_full <- as.integer(a[[2]])
D <- "/home/sina/projects/synchromesh/data"
mj <- jsonlite::fromJSON(file.path(D,ds,"subset_genes",paste0(ds,"_subset_genes_metrics.json")))
stored <- mj$sctransform$cov_gene

M <- t(as(readMM(file.path(D,ds,"subset_genes","raw.mtx.gz")),"CsparseMatrix"))  # genes x cells
rownames(M) <- paste0("g",seq_len(nrow(M))); colnames(M) <- paste0("c",seq_len(ncol(M)))
nc <- ncol(M)
mu  <- Matrix::rowMeans(M); det <- Matrix::rowSums(M>0)/nc
mask <- mu>0.05 & det>=0.01
cat(sprintf("%s: %d genes x %d cells | filtered genes (mu>0.05 & det>=1%%): %d | stored cov_gene=%.3f\n",
            ds, nrow(M), nc, sum(mask), stored))

rowvar <- function(Y){ m<-rowMeans(Y); rowMeans(Y*Y)-m*m }
CVq <- function(v){ v<-v[is.finite(v)&v>0]; sd(v)/mean(v) }
vrun <- function(X) vst(X, vst.flavor="v2", method="glmGamPoi_offset", verbosity=0)$y

kc <- Matrix::colSums(M)>0
if (do_full==1) {
  kg <- Matrix::rowSums(M)>0
  y  <- vrun(M[kg,kc,drop=FALSE])
  vg <- rowvar(y)
  cv_all <- CVq(vg)
  mfilt <- mask[rownames(y)]                       # filtered genes among kept
  cv_full_filt <- CVq(vg[mfilt])
  cat(sprintf("  full-fit          all-genes CV(q)= %.3f   (stored %.3f)\n", cv_all, stored))
  cat(sprintf("  full-fit-filtered          CV(q)= %.3f\n", cv_full_filt))
}
# filtered-fit: re-run vst on the filtered-gene matrix only
gf <- which(mask & Matrix::rowSums(M)>0)
yf <- vrun(M[gf,kc,drop=FALSE])
cat(sprintf("  filtered-fit               CV(q)= %.3f   (genes fit=%d)\n", CVq(rowvar(yf)), nrow(yf)))
