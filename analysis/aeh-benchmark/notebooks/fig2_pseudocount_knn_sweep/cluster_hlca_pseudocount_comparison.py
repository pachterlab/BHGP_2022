#!/usr/bin/env python3
"""Compare Leiden clusterings from PFlog at the estimated and CP10k pseudocounts."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

import scclr


HERE = Path(__file__).resolve().parent
AEH_BENCHMARK_DIR = HERE.parents[1]
OUT_DIR = AEH_BENCHMARK_DIR / "output"
DATA = AEH_BENCHMARK_DIR / "notebooks" / "fig4_lung_depth_geometry" / "data" / "hlca_lung_seqwell_10xv3_balanced_raw.h5ad"
OUT = OUT_DIR / "hlca_lung_balanced_pflog_cp10k_leiden_similarity.tsv"
LABELS = OUT_DIR / "hlca_lung_balanced_pflog_cp10k_leiden_labels.tsv"
PCA_DIM = 50
SEED = 1
N_GENES = 3000
MIN_DETECTION = 0.01


def as_csr(X) -> sp.csr_matrix:
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


def filter_cells_genes(X: sp.csr_matrix) -> tuple[sp.csr_matrix, np.ndarray]:
    cell_keep = np.asarray(X.sum(axis=1)).ravel() > 0
    X = X[cell_keep].tocsr()
    detected = np.asarray((X > 0).mean(axis=0)).ravel()
    gene_keep = detected > MIN_DETECTION
    X = X[:, gene_keep].tocsr()
    score = np.asarray(X.mean(axis=0)).ravel()
    idx = np.argsort(score)[-min(N_GENES, X.shape[1]) :]
    idx = np.sort(idx)
    return X[:, idx].tocsr(), cell_keep


def dense_pca(X_cells_genes: np.ndarray, n_components: int) -> np.ndarray:
    k = min(n_components, X_cells_genes.shape[0] - 1, X_cells_genes.shape[1] - 1)
    return PCA(n_components=int(k), svd_solver="full", random_state=SEED).fit_transform(X_cells_genes)


def pflog_scores(counts: sp.csr_matrix, pseudocount: float, n_components: int) -> np.ndarray:
    depths = np.asarray(counts.sum(axis=1)).ravel()
    sbar = float(depths.mean())
    pf_counts = counts.multiply(sbar / depths[:, None]).toarray()
    Z = np.log(pf_counts + pseudocount)
    Z -= Z.mean(axis=1, keepdims=True)
    return dense_pca(Z, n_components)


def leiden_from_scores(scores: np.ndarray, *, n_neighbors: int) -> np.ndarray:
    a = ad.AnnData(np.zeros((scores.shape[0], 1), dtype=np.float32))
    a.obsm["X_pca"] = scores
    sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X_pca", random_state=SEED)
    sc.tl.leiden(a, random_state=SEED)
    return a.obs["leiden"].astype(str).to_numpy()


def counts_by_label(labels: np.ndarray) -> dict[str, int]:
    return {str(k): int(v) for k, v in pd.Series(labels).value_counts().sort_index().items()}


def main() -> None:
    cache_dir = Path(tempfile.gettempdir()) / "pflog_scanpy_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    adata = ad.read_h5ad(DATA)
    counts, cell_keep = filter_cells_genes(as_csr(adata.X).astype(np.float64))
    obs = adata.obs.iloc[np.flatnonzero(cell_keep)].copy()
    depths = np.asarray(counts.sum(axis=1)).ravel()
    sbar = float(depths.mean())

    pf_res = scclr.normalize_pca(counts, n_components=1, target="auto", seed=SEED, tol=1e-6)
    alpha = float(pf_res.alpha)
    y0_opt = 1.0 / (4.0 * alpha)
    y0_cp10k = sbar / 10_000.0

    opt_scores = pflog_scores(counts, y0_opt, PCA_DIM)
    cp10k_scores = pflog_scores(counts, y0_cp10k, PCA_DIM)

    rows = []
    label_table = pd.DataFrame(
        {
            "cell": obs.index.astype(str),
            "assay": obs["assay"].astype(str).to_numpy() if "assay" in obs else "",
            "cell_type": obs["cell_type"].astype(str).to_numpy() if "cell_type" in obs else "",
        }
    )
    for n_neighbors in (15, 50):
        opt = leiden_from_scores(opt_scores, n_neighbors=n_neighbors)
        cp10k = leiden_from_scores(cp10k_scores, n_neighbors=n_neighbors)
        label_table[f"PFlog_optimum_leiden_k{n_neighbors}"] = opt
        label_table[f"PFlog_CP10k_leiden_k{n_neighbors}"] = cp10k
        rows.append(
            {
                "dataset": "HLCA lung balanced subset",
                "n_cells": int(counts.shape[0]),
                "n_genes": int(counts.shape[1]),
                "pca_dim": PCA_DIM,
                "neighbors": n_neighbors,
                "alpha": alpha,
                "PFlog_optimum_pseudocount": y0_opt,
                "CP10k_pseudocount": y0_cp10k,
                "adjusted_rand_index": adjusted_rand_score(opt, cp10k),
                "normalized_mutual_information": normalized_mutual_info_score(opt, cp10k),
                "PFlog_optimum_n_clusters": int(len(set(opt))),
                "CP10k_n_clusters": int(len(set(cp10k))),
                "PFlog_optimum_cluster_sizes": json.dumps(counts_by_label(opt), sort_keys=True),
                "CP10k_cluster_sizes": json.dumps(counts_by_label(cp10k), sort_keys=True),
            }
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    label_table.to_csv(LABELS, sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"wrote: {OUT}")
    print(f"wrote: {LABELS}")


if __name__ == "__main__":
    main()
