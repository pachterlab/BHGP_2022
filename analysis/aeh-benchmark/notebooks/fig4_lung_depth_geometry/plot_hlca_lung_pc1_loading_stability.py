#!/usr/bin/env python3
"""Generate main Fig. 3: cross-depth PC1 loading stability in lung cells."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "data" / "hlca_lung_seqwell_10xv3_balanced_raw.h5ad"
DEFAULT_FIGURE = HERE / "figures" / "lung_10x_vs_seqwell_depth_stability_pc1_loadings.pdf"
DEFAULT_TABLE = HERE / "tables" / "lung_10x_vs_seqwell_depth_stability_summary.tsv"
DEFAULT_LOADINGS = HERE / "tables" / "lung_10x_vs_seqwell_depth_stability_pc1_loadings.tsv"


def import_scclr(extra_path: Path | None = None):
    if extra_path is not None:
        sys.path.insert(0, str(extra_path))
    import scclr

    return scclr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--summary", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--loadings", type=Path, default=DEFAULT_LOADINGS)
    parser.add_argument("--scclr-path", type=Path, default=None, help="Optional path to a local scclr/python directory.")
    parser.add_argument("--n-genes", type=int, default=1000)
    parser.add_argument("--min-detection", type=float, default=0.01)
    return parser.parse_args()


def as_csr(X) -> sp.csr_matrix:
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


def log_cp10k(X: sp.csr_matrix) -> np.ndarray:
    depth = np.asarray(X.sum(axis=1)).ravel()
    scale = np.divide(10000.0, depth, out=np.zeros_like(depth, dtype=float), where=depth > 0)
    Y = X.multiply(scale[:, None]).tocsr()
    Y.data = np.log1p(Y.data)
    return Y.toarray()


def logpf_with_k(X: sp.csr_matrix, k: float) -> np.ndarray:
    depth = np.asarray(X.sum(axis=1)).ravel()
    scale = np.divide(float(k), depth, out=np.zeros_like(depth, dtype=float), where=depth > 0)
    Y = X.multiply(scale[:, None]).tocsr()
    Y.data = np.log1p(Y.data)
    return Y.toarray()


def pca_pc1_loadings_dense(Y: np.ndarray) -> np.ndarray:
    Y = Y - Y.mean(axis=0, keepdims=True)
    model = PCA(n_components=1, svd_solver="full")
    model.fit(Y)
    return model.components_[0].copy()


def pca_pc1_loadings_pflog(scclr, X: sp.csr_matrix) -> tuple[np.ndarray, dict[str, float]]:
    res = scclr.normalize_pca(X, n_components=1, target="auto", seed=1)
    return np.asarray(res.components)[0].copy(), {"alpha": float(res.alpha), "k": float(res.k)}


def align(reference: np.ndarray, candidate: np.ndarray) -> tuple[np.ndarray, float]:
    r = float(np.corrcoef(reference, candidate)[0, 1])
    if r < 0:
        return -candidate, -r
    return candidate, r


def top_common_genes(
    X_low: sp.csr_matrix,
    X_high: sp.csr_matrix,
    genes: np.ndarray,
    n_top: int,
    min_detection: float,
) -> tuple[sp.csr_matrix, sp.csr_matrix, np.ndarray]:
    keep = (np.asarray((X_low > 0).mean(axis=0)).ravel() > min_detection) & (
        np.asarray((X_high > 0).mean(axis=0)).ravel() > min_detection
    )
    X_low = X_low[:, keep].tocsr()
    X_high = X_high[:, keep].tocsr()
    genes = genes[keep]
    score = np.asarray(X_low.mean(axis=0)).ravel() + np.asarray(X_high.mean(axis=0)).ravel()
    idx = np.argsort(score)[-min(n_top, len(score)) :]
    idx = idx[np.argsort(genes[idx])]
    return X_low[:, idx].tocsr(), X_high[:, idx].tocsr(), genes[idx]


def main() -> None:
    args = parse_args()
    scclr = import_scclr(args.scclr_path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.loadings.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    genes = adata.var["feature_name"].astype(str).to_numpy() if "feature_name" in adata.var else adata.var_names.astype(str).to_numpy()
    X = as_csr(adata.X).astype(np.float64).tocsr()
    obs = adata.obs.copy()
    assay = obs["assay"].astype(str).to_numpy()
    low_mask = assay == "Seq-Well"
    high_mask = assay == "10x 3' v3"
    X_low_full = X[low_mask].tocsr()
    X_high_full = X[high_mask].tocsr()
    low_depth = np.asarray(X_low_full.sum(axis=1)).ravel()
    high_depth = np.asarray(X_high_full.sum(axis=1)).ravel()
    X_low, X_high, genes = top_common_genes(
        X_low_full,
        X_high_full,
        genes,
        n_top=args.n_genes,
        min_detection=args.min_detection,
    )

    low_pf, low_norm = pca_pc1_loadings_pflog(scclr, X_low)
    high_pf, high_norm = pca_pc1_loadings_pflog(scclr, X_high)
    comparisons = {
        "logCP10K": (pca_pc1_loadings_dense(log_cp10k(X_low)), pca_pc1_loadings_dense(log_cp10k(X_high)), {}, {}),
        "logPF, estimated K": (
            pca_pc1_loadings_dense(logpf_with_k(X_low, low_norm["k"])),
            pca_pc1_loadings_dense(logpf_with_k(X_high, high_norm["k"])),
            low_norm,
            high_norm,
        ),
        "PFlog": (low_pf, high_pf, low_norm, high_norm),
    }

    loadings = pd.DataFrame({"gene": genes})
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7), constrained_layout=True)
    for panel, (ax, (method, (l_low, l_high, low_method_norm, high_method_norm))) in zip("abc", zip(axes, comparisons.items())):
        l_high, r = align(l_low, l_high)
        r2 = r**2
        loadings[f"{method}_Seq-Well_PC1_loading"] = l_low
        loadings[f"{method}_10x_v3_PC1_loading"] = l_high
        rows.append(
            {
                "method": method,
                "pc1_loading_correlation_abs_sign_aligned": r,
                "pc1_loading_r2_abs_sign_aligned": r2,
                "n_genes": len(genes),
                "seqwell_n_cells": int(X_low_full.shape[0]),
                "seqwell_mean_depth": float(low_depth.mean()),
                "seqwell_median_depth": float(np.median(low_depth)),
                "seqwell_alpha": low_method_norm.get("alpha", np.nan),
                "seqwell_K": low_method_norm.get("k", np.nan),
                "tenx_v3_n_cells": int(X_high_full.shape[0]),
                "tenx_v3_mean_depth": float(high_depth.mean()),
                "tenx_v3_median_depth": float(np.median(high_depth)),
                "tenx_v3_alpha": high_method_norm.get("alpha", np.nan),
                "tenx_v3_K": high_method_norm.get("k", np.nan),
            }
        )
        lim = np.nanmax(np.abs(np.r_[l_low, l_high]))
        ax.scatter(l_low, l_high, s=7, alpha=0.35, linewidths=0, color="#2f5d7c")
        ax.plot([-lim, lim], [-lim, lim], color="#9a3324", lw=1.1)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_title(f"{method}: $R^2$ = {r2:.2f}")
        ax.set_xlabel("lung Seq-Well PC1 loading")
        ax.set_ylabel("lung 10x v3 PC1 loading")
        ax.text(-0.12, 1.08, panel, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")

    summary = pd.DataFrame(rows)
    loadings.to_csv(args.loadings, sep="\t", index=False)
    summary.to_csv(args.summary, sep="\t", index=False)
    fig.savefig(args.output)
    fig.savefig(args.output.with_suffix(".png"), dpi=300)
    print(summary.to_string(index=False))
    print(f"wrote: {args.output}")


if __name__ == "__main__":
    main()
