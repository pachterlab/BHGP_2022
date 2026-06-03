#!/usr/bin/env python3
"""Compare sctransform implementations (pysctransform 0.1.1 vs R v2/glmGamPoi)
on the manuscript's metrics for one dataset, to decide whether claims that
depend on sctransform need to be regenerated.

Computes, for each implementation's residual matrix:
  - global (summary-boxplot) metrics: cov_gene, r2_depth, r_mono
  - celltype-bar metrics on a chosen cell type: pc1_entropy_frac, fp_de_genes,
    mean_abs_spearman

Usage:
  python compare_sct_impl.py <dataset_dir> <pysct_csv_gz> <rsct_csv_gz> <celltype>
"""
import gzip
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse, stats
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA


def load_csv(path):
    with gzip.open(path, "rt") as f:
        return np.loadtxt(f, delimiter=",", dtype=np.float32)


def load_any(path):
    """CSV(.gz) for pysctransform, or float32 binary (+ .shape sidecar) for R."""
    if path.endswith(".bin"):
        n_cells, n_genes = map(int, open(path + ".shape").read().strip().split(","))
        return np.fromfile(path, dtype=np.float32).reshape(n_cells, n_genes)
    return load_csv(path)


# ---- summary metrics (from metrics_methods_dense.py) ----
def cov_gene(mtx):
    gvar = np.var(mtx, axis=0)
    return float(np.sqrt(np.var(gvar)) / np.mean(gvar))


def r2_depth(mtx, raw):
    x = np.asarray(raw.sum(1)).ravel().astype(float)
    y = mtx.sum(1).ravel().astype(float)
    minx, maxx = x.min(), x.max()
    miny, maxy = y.min(), y.max()
    maxy = maxy - miny
    if maxy == 0:
        return 0.0
    xx = (x - minx) / maxx
    yy = (y - miny) / maxy
    reg = LinearRegression().fit(xx.reshape(-1, 1), yy)
    return float(reg.score(xx.reshape(-1, 1), yy))


def r_mono(mtx, raw):
    rv = np.empty(mtx.shape[0])
    for i in range(mtx.shape[0]):
        rv[i], _ = stats.spearmanr(mtx[i], raw[i].toarray().ravel())
    return float(np.nanmean(rv))


# ---- celltype metrics (from metrics_celltype.py) ----
def pc1_entropy_frac(X):
    pca = PCA(n_components=1, svd_solver="full").fit(X)
    w = np.abs(pca.components_[0]); s = w.sum()
    if s <= 0:
        return float("nan")
    p = w / s
    return float(stats.entropy(p) / np.log(len(p)))


def fp_de_genes(X_full, X, raw_depth, n_top=500, alpha=0.01, nan_cutoff=0.1):
    order = np.argsort(raw_depth)
    n_top = min(n_top, len(order) // 2)
    bot, top = order[:n_top], order[-n_top:]
    raw_top = X_full[top]
    keep = (raw_top > 0).sum(0) > nan_cutoff * raw_top.shape[0]
    if keep.sum() == 0:
        return 0
    Xt, Xb = X[top][:, keep], X[bot][:, keep]
    res = stats.ttest_ind(Xt, Xb, axis=0, equal_var=False, nan_policy="propagate")
    p = np.where(np.isnan(res.pvalue), 1.0, res.pvalue) / 2.0
    return int((np.minimum(1.0, p * keep.sum()) < alpha).sum())


def mean_abs_spearman(X, n_sub=500, seed=0):
    rng = np.random.default_rng(seed)
    if X.shape[0] > n_sub:
        X = X[rng.choice(X.shape[0], n_sub, replace=False)]
    X = X.astype(np.float64, copy=True)
    flat = X.ravel(); nz = flat[flat != 0]
    if nz.size:
        u = np.unique(np.concatenate(([0.0], nz))); d = np.diff(u); d = d[d > 0]
        md = float(d.min()) if d.size else 1.0
    else:
        md = 1.0
    X = X + rng.uniform(-md / 4, md / 4, size=X.shape[1])
    rho = stats.spearmanr(X, axis=1).correlation
    iu = np.triu_indices_from(rho, k=1)
    return float(np.mean(np.abs(rho[iu])))


def report(name, X, raw_full, raw_depth_ct, idx, X_full_ct):
    print(f"\n=== {name} ===")
    print(f"  cov_gene  = {cov_gene(X):.4f}")
    print(f"  r2_depth  = {r2_depth(X, raw_full):.4f}")
    print(f"  r_mono    = {r_mono(X, raw_full):.4f}")
    Xct = X[idx]
    print(f"  [celltype] pc1_entropy_frac = {pc1_entropy_frac(Xct):.4f}")
    print(f"  [celltype] fp_de_genes      = {fp_de_genes(X_full_ct, Xct, raw_depth_ct)}")
    print(f"  [celltype] mean_abs_spearman= {mean_abs_spearman(Xct):.4f}")


def main(dataset_dir, pysct_fn, rsct_fn, celltype):
    raw = mmread(f"{dataset_dir}/subset_genes/raw.mtx.gz").tocsr().astype(np.float32)
    meta = pd.read_csv(f"{dataset_dir}/metadata_barcodes.txt.gz")
    ct = meta[meta.columns[-1]].values
    idx = np.where(ct == celltype)[0]
    print(f"raw: {raw.shape}; celltype {celltype!r}: {idx.size} cells")
    raw_ct = np.asarray(raw[idx].todense())
    raw_depth_ct = raw_ct.sum(1)

    for name, fn in [("pysctransform 0.1.1", pysct_fn), ("R sctransform v2", rsct_fn)]:
        X = load_any(fn)
        report(name, X, raw, raw_depth_ct, idx, raw_ct)


if __name__ == "__main__":
    main(*sys.argv[1:5])
