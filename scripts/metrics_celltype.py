#!/usr/bin/env python3
"""Compute three cell-type-level metrics for each normalization method on a
specified cell type within a dataset.

Metrics (per method):
  (a) pc1_entropy_frac : entropy of |l1-normalized PC1 loadings| / log(n_genes)
  (b) fp_de_genes      : # false-positive DE genes between top-500 and bottom-500
                         cells by raw depth (Welch t-test on per-gene values,
                         BH-corrected p<0.01)
  (c) mean_abs_spearman: mean of |pairwise Spearman r| between cells (with
                         tie-breaking jitter; sub-sampled to N cells for speed)

Usage:
    python metrics_celltype.py <dataset_dir> <celltype> <out_json> [--n_sub 500]

<dataset_dir> must contain:
    matrix.mtx.gz, barcodes.txt.gz, metadata_barcodes.txt.gz (cells × genes raw)
    subset_genes/{raw,pf,sqrt,log,cp10k_log,cpm_log,pf_log,pf_log_pf}.mtx.gz
    subset_genes/{sctransform,cp10k_log_scale}.csv.gz
The metadata file's last column must be 'celltype' and the first column the barcode.
"""
import argparse
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread
from sklearn.decomposition import PCA

# Method label → input file under subset_genes/, and whether it's dense.
METHODS = [
    ("raw",                  "raw.mtx.gz",            False),
    ("PF",                   "pf.mtx.gz",             False),
    ("sqrt",                 "sqrt.mtx.gz",           False),
    ("log1p",                "log.mtx.gz",            False),
    ("log1pCP10k",           "cp10k_log.mtx.gz",      False),
    ("log1pCPM",             "cpm_log.mtx.gz",        False),
    ("scalelog1pCP10k",      "cp10k_log_scale.csv.gz", True),
    ("sctransform",          "sctransform.csv.gz",    True),
    ("log1pPF",              "pf_log.mtx.gz",         False),
    ("PFlogPF (shift. CLR)", "pf_log_pf.mtx.gz",      False),
]


def load_matrix(path, dense):
    """Return a (cells × genes) dense float32 numpy array."""
    if dense:
        with gzip.open(path, "rt") as f:
            arr = np.loadtxt(f, delimiter=",", dtype=np.float32)
    else:
        m = mmread(path)
        arr = np.asarray(m.todense(), dtype=np.float32) if sparse.issparse(m) else np.asarray(m, dtype=np.float32)
    return arr


def pc1_entropy_frac(X):
    """Entropy of |l1-normalized PC1 loadings| / log(n_genes)."""
    pca = PCA(n_components=1, svd_solver="full")
    pca.fit(X)
    w = np.abs(pca.components_[0])
    s = w.sum()
    if s <= 0:
        return float("nan")
    p = w / s
    ent = stats.entropy(p)
    max_ent = np.log(len(p))
    return float(ent / max_ent) if max_ent > 0 else float("nan")


def fp_de_genes(X_full, X, raw_depth, n_top=500, alpha=0.01, nan_cutoff=0.1):
    """# DE genes with Bonferroni-corrected p < alpha between top-n_top and
    bottom-n_top cells by raw depth. Matches Booeshaghi et al. 2021 dexpress:
    Welch t-test, gene filter requires >nan_cutoff*ncells nonzero entries in
    the target group, p_raw halved (one-sided convention), Bonferroni-corrected.

    X_full is the *raw* counts matrix (used for the >nan_cutoff expression
    filter so the gene-passing set is consistent across methods); X is the
    method's normalized matrix (used for the actual t-test)."""
    order = np.argsort(raw_depth)
    n_top = min(n_top, len(order) // 2)
    bot = order[:n_top]
    top = order[-n_top:]
    # Gene filter on raw counts: keep genes expressed in > nan_cutoff fraction
    # of the target (top) group.
    raw_top = X_full[top]
    keep = ((raw_top > 0).sum(axis=0) > nan_cutoff * raw_top.shape[0])
    if keep.sum() == 0:
        return 0
    Xt, Xb = X[top][:, keep], X[bot][:, keep]
    res = stats.ttest_ind(Xt, Xb, axis=0, equal_var=False, nan_policy="propagate")
    p = np.asarray(res.pvalue, dtype=float)
    p = np.where(np.isnan(p), 1.0, p) / 2.0   # one-sided convention
    # Bonferroni
    p_adj = np.minimum(1.0, p * keep.sum())
    return int((p_adj < alpha).sum())


def mean_abs_spearman(X, n_sub=500, seed=0):
    """Mean of |pairwise Spearman r| between cells with per-gene tie-break
    offset, sub-sampled to n_sub cells if larger.

    Matches the original angelidis_mono.ipynb pairwise_spearman: the offset
    is a single random shift per gene (constant across cells), so ties
    between cells at the same gene break consistently — the rank-order of
    'all zeros' positions becomes the same in every cell, contributing
    correctly to the cross-cell correlation."""
    rng = np.random.default_rng(seed)
    if X.shape[0] > n_sub:
        idx = rng.choice(X.shape[0], n_sub, replace=False)
        X = X[idx]
    X = X.astype(np.float64, copy=True)
    # Smallest nonzero difference between sorted unique entries.
    flat = X.ravel()
    nz = flat[flat != 0]
    if nz.size > 0:
        u = np.unique(np.concatenate(([0.0], nz)))
        diffs = np.diff(u)
        diffs = diffs[diffs > 0]
        min_diff = float(diffs.min()) if diffs.size > 0 else 1.0
    else:
        min_diff = 1.0
    # Per-gene constant offset, range [-min_diff/4, min_diff/4].
    offsets = rng.uniform(-min_diff / 4, min_diff / 4, size=X.shape[1])
    X = X + offsets  # broadcasts across rows
    # spearmanr with axis=1 treats each ROW as a variable -> correlation
    # between cells.
    rho = stats.spearmanr(X, axis=1).correlation
    if np.ndim(rho) == 0:
        return float(abs(rho))
    iu = np.triu_indices_from(rho, k=1)
    return float(np.mean(np.abs(rho[iu])))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dataset_dir")
    ap.add_argument("celltype")
    ap.add_argument("out_json")
    ap.add_argument("--n_sub", type=int, default=500,
                    help="cells to subsample for pairwise Spearman (default 500)")
    args = ap.parse_args()

    # Load metadata, find cell-type indices in matrix-column order.
    meta_path = os.path.join(args.dataset_dir, "metadata_barcodes.txt.gz")
    meta = pd.read_csv(meta_path)
    barcode_col = meta.columns[0]
    celltype_col = meta.columns[-1]
    print(f"metadata: {meta.shape[0]} cells, celltype col = {celltype_col!r}", file=sys.stderr)

    # The matrix's cell order is the metadata's row order.
    is_target = (meta[celltype_col] == args.celltype).values
    idx = np.where(is_target)[0]
    if idx.size == 0:
        sys.exit(f"no cells with celltype == {args.celltype!r}")
    print(f"{args.celltype}: {idx.size} cells", file=sys.stderr)

    # Raw depth (per-cell sum on raw counts) — same for all methods, derived
    # from subset_genes/raw.mtx.gz.
    raw_path = os.path.join(args.dataset_dir, "subset_genes", "raw.mtx.gz")
    print("loading raw for depth...", file=sys.stderr)
    raw_full = load_matrix(raw_path, dense=False)
    raw_sub = raw_full[idx]
    raw_depth = raw_sub.sum(axis=1)

    results = {
        "dataset": os.path.basename(os.path.normpath(args.dataset_dir)),
        "celltype": args.celltype,
        "n_cells": int(idx.size),
        "n_genes": int(raw_full.shape[1]),
        "methods": {},
    }

    for label, fname, dense in METHODS:
        path = os.path.join(args.dataset_dir, "subset_genes", fname)
        if not os.path.exists(path):
            print(f"  {label}: skip (missing {fname})", file=sys.stderr)
            continue
        print(f"  {label}: load+compute...", file=sys.stderr)
        X = load_matrix(path, dense=dense)
        X_sub = X[idx]
        del X  # free full matrix
        pc1 = pc1_entropy_frac(X_sub)
        fp = fp_de_genes(raw_sub, X_sub, raw_depth)
        ms = mean_abs_spearman(X_sub, n_sub=args.n_sub)
        results["methods"][label] = {
            "pc1_entropy_frac": pc1,
            "fp_de_genes": fp,
            "mean_abs_spearman": ms,
        }
        print(f"    pc1={pc1:.4f}  fp_de={fp}  |spearman|={ms:.4f}",
              file=sys.stderr)

    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"wrote {args.out_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
