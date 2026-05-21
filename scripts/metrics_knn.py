#!/usr/bin/env python3
"""k-NN consistency metric (Ahlmann-Eltze & Huber, Nat Methods 2023).

For each cell, the metric is the overlap of its k nearest neighbors when the
gene set is randomly split in half and the transformation is applied to each
half independently. A good transformation gives consistent k-NN graphs across
random gene subsets.

  overlap_per_cell(i) = |knn_A(i) intersect knn_B(i)| / k
  k-NN_consistency  = mean over cells

For each method we compute the metric at fixed defaults:
  - k = 50 neighbors (AE&H default)
  - PCA dimension = 50 (AE&H tested 5-1000; 50 is their main figure)
  - seed = 42 (fixed for reproducibility)

Usage:
    python metrics_knn.py <raw.mtx.gz> <method_name> [<output.json>]

If output.json is omitted, prints to stdout.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.io import mmread
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from norm_sparse import (  # noqa: E402
    norm_raw, norm_pf, norm_log, norm_pf_log, norm_pf_log_pf,
    norm_cpm_log, norm_cp10k_log, norm_sqrt, norm_clr,
)

METHODS = {
    "raw":        norm_raw,
    "pf":         norm_pf,
    "log":        norm_log,
    "pf_log":     norm_pf_log,
    "pf_log_pf":  norm_pf_log_pf,
    "cpm_log":    norm_cpm_log,
    "cp10k_log":  norm_cp10k_log,
    "sqrt":       norm_sqrt,
    "clr":        norm_clr,
}


def _resolve_method(name, c=None):
    """Allow 'clr@<c>' syntax to pass a pseudocount, e.g. 'clr@0.4'."""
    if "@" in name:
        base, c_str = name.split("@", 1)
        c_val = float(c_str)
    else:
        base = name
        c_val = c
    fn = METHODS[base]
    if base == "clr" and c_val is not None:
        return lambda m, _fn=fn, _c=c_val: _fn(m, c=_c), base
    return fn, base


def _ensure_dense(m):
    return m.toarray() if sp.issparse(m) else m


def knn_overlap_for_method(raw, method_name, k=50, n_pca=50, seed=42):
    """Compute split-half k-NN overlap for one transformation.

    method_name accepts 'clr@<c>' syntax to set pseudocount, e.g. 'clr@0.4032'.
    """
    n_cells, n_genes = raw.shape
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_genes)
    half_a = perm[: n_genes // 2]
    half_b = perm[n_genes // 2 :]

    fn, _ = _resolve_method(method_name)

    t0 = time.time()
    out_a = fn(raw[:, half_a])
    out_b = fn(raw[:, half_b])
    t_transform = time.time() - t0

    t0 = time.time()
    out_a = _ensure_dense(out_a)
    out_b = _ensure_dense(out_b)
    t_densify = time.time() - t0

    n_comp = min(n_pca, n_cells - 1, out_a.shape[1] - 1, out_b.shape[1] - 1)
    if n_comp <= 0:
        return float("nan")

    t0 = time.time()
    pca_a = PCA(n_components=n_comp, random_state=seed).fit_transform(out_a)
    pca_b = PCA(n_components=n_comp, random_state=seed).fit_transform(out_b)
    t_pca = time.time() - t0

    k_eff = min(k, n_cells - 1)
    t0 = time.time()
    knn_a = NearestNeighbors(n_neighbors=k_eff + 1).fit(pca_a).kneighbors(
        pca_a, return_distance=False
    )
    knn_b = NearestNeighbors(n_neighbors=k_eff + 1).fit(pca_b).kneighbors(
        pca_b, return_distance=False
    )
    t_knn = time.time() - t0

    t0 = time.time()
    overlaps = np.empty(n_cells)
    for i in range(n_cells):
        overlaps[i] = len(set(knn_a[i, 1:]).intersection(knn_b[i, 1:])) / k_eff
    t_overlap = time.time() - t0

    return {
        "k": k_eff,
        "n_pca": n_comp,
        "seed": seed,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "mean_overlap": float(overlaps.mean()),
        "median_overlap": float(np.median(overlaps)),
        "timings": {
            "transform_s": round(t_transform, 2),
            "densify_s": round(t_densify, 2),
            "pca_s": round(t_pca, 2),
            "knn_s": round(t_knn, 2),
            "overlap_s": round(t_overlap, 2),
        },
    }


def main(raw_path, method_name, out_path=None, k=50, n_pca=50, seed=42):
    print(f"loading {raw_path}", flush=True)
    raw = mmread(raw_path).tocsr()
    print(f"  shape: {raw.shape}", flush=True)

    print(f"k-NN overlap for {method_name}  (k={k}, n_pca={n_pca}, seed={seed})", flush=True)
    result = knn_overlap_for_method(raw, method_name, k=k, n_pca=n_pca, seed=seed)
    print(json.dumps(result, indent=2))

    if out_path:
        if os.path.exists(out_path):
            with open(out_path) as f:
                doc = json.load(f)
        else:
            doc = {}
        doc.setdefault("knn_overlap", {})[method_name] = result
        with open(out_path, "w") as f:
            json.dump(doc, f, indent=2)
        print(f"wrote {out_path}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    raw_path = sys.argv[1]
    method = sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else None
    main(raw_path, method, out)
