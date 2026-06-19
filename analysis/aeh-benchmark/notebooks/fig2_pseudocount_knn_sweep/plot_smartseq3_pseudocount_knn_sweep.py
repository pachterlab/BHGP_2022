#!/usr/bin/env python3
"""Sweep PFlog pseudocounts and compare kNN graphs to PFlog on a Fig. 2c dataset."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

import scclr


matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
AEH_BENCHMARK_DIR = HERE.parents[1]
REPO = AEH_BENCHMARK_DIR.parent
DATA = AEH_BENCHMARK_DIR / "benchmark" / "output" / "clr_local" / "data" / "downsampling"
OUT_DIR = AEH_BENCHMARK_DIR / "output"
OUT = OUT_DIR / "smartseq3_fibroblasts_pflog_pseudocount_knn_overlap.pdf"
DATASET = "smartSeq3_fibroblasts"
PCA_DIM = 10
KNN = 50
SEED = 1


def read_table_counts(path: Path) -> sp.csr_matrix:
    df = pd.read_csv(path, sep="\t", index_col=0)
    return sp.csr_matrix(df.to_numpy(dtype=np.float64))


def filter_gene_cell(M_gene_cell: sp.spmatrix) -> sp.csr_matrix:
    M = M_gene_cell.tocsr()
    gene_keep = np.asarray(M.sum(axis=1)).ravel() > 0
    cell_keep = np.asarray(M.sum(axis=0)).ravel() > 0
    return M[gene_keep][:, cell_keep].tocsr()


def knn(emb: np.ndarray, k: int) -> np.ndarray:
    kk = min(k + 1, emb.shape[0])
    return NearestNeighbors(n_neighbors=kk).fit(emb).kneighbors(return_distance=False)[:, 1:]


def mean_overlap(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean([len(set(x).intersection(y)) for x, y in zip(a, b)]))


def dense_pca(X_cells_genes: np.ndarray, n_components: int) -> np.ndarray:
    k = min(n_components, X_cells_genes.shape[0] - 1, X_cells_genes.shape[1] - 1)
    return PCA(n_components=int(k), svd_solver="full").fit_transform(X_cells_genes)


def pflog_transform(pf_counts: np.ndarray, pseudocount: float) -> np.ndarray:
    Z = np.log(pf_counts + pseudocount)
    return Z - Z.mean(axis=1, keepdims=True)


def main() -> None:
    cache_dir = Path(tempfile.gettempdir()) / "pflog_matplotlib_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    counts_gene_cell = filter_gene_cell(read_table_counts(DATA / "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt"))
    counts = counts_gene_cell.T.tocsr().astype(np.float64)
    depths = np.asarray(counts.sum(axis=1)).ravel()
    sbar = float(depths.mean())

    pf_res = scclr.normalize_pca(counts, n_components=PCA_DIM, target="auto", seed=SEED, tol=1e-6)
    pf_scores = np.asarray(pf_res.scores)
    ref_knn = knn(pf_scores, KNN)
    alpha = float(pf_res.alpha)
    k_auto = float(pf_res.k)
    y0_opt = 1.0 / (4.0 * alpha)

    reference = {
        "CPM": sbar / 1_000_000.0,
        "log1pPF": 1.0,
        "PFlog optimum": y0_opt,
        "CP10k": sbar / 10_000.0,
    }
    lo = min(reference.values()) / 8.0
    hi = max(reference.values()) * 8.0
    pseudocounts = np.unique(np.r_[np.logspace(np.log10(lo), np.log10(hi), 95), list(reference.values())])

    pf_counts = counts.multiply(sbar / depths[:, None]).toarray()
    rows = []
    for y0 in pseudocounts:
        emb = dense_pca(pflog_transform(pf_counts, float(y0)), PCA_DIM)
        rows.append(
            {
                "dataset": DATASET,
                "pseudocount": float(y0),
                "shared_neighbors_with_pflog": mean_overlap(ref_knn, knn(emb, KNN)),
                "k": KNN,
                "pca_dim": PCA_DIM,
            }
        )
    sweep = pd.DataFrame(rows).sort_values("pseudocount")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sweep.to_csv(OUT.with_suffix(".tsv"), sep="\t", index=False)

    matplotlib.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
        }
    )
    fig, ax = plt.subplots(figsize=(4.8, 3.2), constrained_layout=True)
    ax.plot(
        sweep["pseudocount"],
        sweep["shared_neighbors_with_pflog"],
        color="#4d4d4d",
        lw=1.9,
        label="PFlog",
    )

    marker_styles = {
        "PFlog optimum": dict(color="#e41a1c", linestyle="-", lw=1.0),
        "CPM": dict(color="#2c7fb8", linestyle="--", lw=0.9),
        "log1pPF": dict(color="#238b45", linestyle="--", lw=0.9),
        "CP10k": dict(color="#756bb1", linestyle="--", lw=0.9),
    }
    label_y = {
        "CPM": 15.0,
        "log1pPF": 19.0,
        "PFlog optimum": 23.0,
        "CP10k": 17.0,
    }
    for name, value in reference.items():
        style = marker_styles[name]
        ax.axvline(value, **style)
        y = np.interp(np.log10(value), np.log10(sweep["pseudocount"]), sweep["shared_neighbors_with_pflog"])
        ax.scatter([value], [y], s=18, color=style["color"], zorder=4)
        ax.text(
            value,
            label_y[name],
            f"{name}\n{value:.3g}",
            rotation=90,
            va="bottom",
            ha="right",
            color=style["color"],
            fontsize=6.8,
        )

    ax.set_xscale("log")
    ax.set_ylim(0, KNN)
    ax.set_xlabel("count-scale pseudocount")
    ax.set_ylabel(f"shared neighbors with PFlog kNN (k={KNN})")
    ax.legend(frameon=False, loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.7, length=3)

    fig.savefig(OUT, bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    summary = {
        "dataset": DATASET,
        "n_cells": int(counts.shape[0]),
        "n_genes": int(counts.shape[1]),
        "pca_dim": PCA_DIM,
        "k": KNN,
        "alpha": alpha,
        "scclr_K": k_auto,
        "mean_depth": sbar,
        "reference_pseudocounts": reference,
        "outputs": {
            "pdf": str(OUT),
            "png": str(OUT.with_suffix(".png")),
            "tsv": str(OUT.with_suffix(".tsv")),
            "json": str(OUT.with_suffix(".json")),
        },
    }
    OUT.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
