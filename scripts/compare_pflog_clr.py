#!/usr/bin/env python3
"""Compare legacy PFlog-ish artifacts against current scclr PFlog.

The legacy ``pf_log_pf.mtx.gz`` artifact is kept only for historical comparison.
Current PFlog is computed on raw counts with scclr as
``center(log1p(4*alpha*x))``.

Reports the summary-boxplot metrics (cov_gene, r2_depth, r_mono) and the
celltype-bar metrics (pc1_entropy_frac, fp_de_genes, mean_abs_spearman).

Usage: python compare_pflog_clr.py <dataset_dir> <celltype>
"""
import sys
import numpy as np, pandas as pd
from scipy.io import mmread
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import compare_sct_impl as C   # reuse identical metric functions
from metrics_matrix import compute_overdispersion
from scclr_pflog import normalize_pflog, to_dense


def main(dsdir, celltype):
    raw = mmread(f"{dsdir}/subset_genes/raw.mtx.gz").tocsr().astype(np.float32)
    mult = np.asarray(mmread(f"{dsdir}/subset_genes/pf_log_pf.mtx.gz").todense(), np.float32)
    alpha = float(compute_overdispersion(raw))
    sclr = to_dense(normalize_pflog(raw, alpha=alpha))
    meta = pd.read_csv(f"{dsdir}/metadata_barcodes.txt.gz")
    idx = np.where(meta[meta.columns[-1]].values == celltype)[0]
    raw_ct = np.asarray(raw[idx].todense()); depth = raw_ct.sum(1)
    print(f"raw {raw.shape} | legacy {mult.shape} | scclr {sclr.shape} | {celltype}: {idx.size} cells\n")

    for name, X in [("legacy pf_log_pf artifact", mult), ("scclr PFlog", sclr)]:
        print(f"=== {name} ===")
        print(f"  row-sum |max|         = {np.abs(X.sum(1)).max():.3f}   (0 => centered/CLR)")
        print(f"  [summary] cov_gene    = {C.cov_gene(X):.4f}")
        print(f"  [summary] r2_depth    = {C.r2_depth(X, raw):.4f}")
        print(f"  [summary] r_mono      = {C.r_mono(X, raw):.4f}")
        Xc = X[idx]
        print(f"  [celltype] pc1_entropy_frac = {C.pc1_entropy_frac(Xc):.4f}")
        print(f"  [celltype] fp_de_genes      = {C.fp_de_genes(raw_ct, Xc, depth)}")
        print(f"  [celltype] mean_abs_spearman= {C.mean_abs_spearman(Xc):.4f}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
