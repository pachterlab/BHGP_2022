#!/usr/bin/env python3
"""Sweep count-scale pseudocounts for the Angelidis pseudobulk PC1/LFC test."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr
from sklearn.decomposition import PCA


matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
AEH_BENCHMARK_DIR = HERE.parents[1]
OUT_DIR = AEH_BENCHMARK_DIR / "output"
DEFAULT_INPUT = HERE / "angelidis_lung_pseudobulk.h5ad"
DEFAULT_DE = HERE / "pc1_loadings_with_edgepython_de.csv"
DEFAULT_OUT = OUT_DIR / "angelidis_pseudobulk_pseudocount_r2_sweep.pdf"


def as_csr(x) -> sp.csr_matrix:
    return x.tocsr() if sp.issparse(x) else sp.csr_matrix(np.asarray(x, dtype=np.float64))


def estimate_alpha(counts: sp.csr_matrix) -> tuple[float, float]:
    depths = np.asarray(counts.sum(axis=1)).ravel()
    sbar = float(depths.mean())
    sf = depths / sbar
    pf = counts.multiply(1.0 / sf[:, None]).toarray()
    mu = pf.mean(axis=0)
    keep = mu > 0
    pf = pf[:, keep]
    mu = mu[keep]
    var = pf.var(axis=0, ddof=1)
    poisson = mu * float(np.mean(1.0 / sf))
    alpha_g = (var - poisson) / (mu * mu)
    alpha_g = alpha_g[np.isfinite(alpha_g) & (alpha_g > 0)]
    if alpha_g.size == 0:
        raise ValueError("could not estimate positive overdispersion")
    return float(np.median(alpha_g)), sbar


def orient_old_positive(scores: np.ndarray, loadings: np.ndarray, ages: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    old = ages == "Old"
    young = ages == "Young"
    if scores[old, 0].mean() < scores[young, 0].mean():
        scores = scores.copy()
        loadings = loadings.copy()
        scores[:, 0] *= -1.0
        loadings[:, 0] *= -1.0
    return scores, loadings


def pc1_loadings_for_pseudocount(
    pf_counts: np.ndarray,
    ages: np.ndarray,
    pseudocount: float,
    clr_center: bool,
) -> np.ndarray:
    transformed = np.log(pf_counts + float(pseudocount))
    if clr_center:
        transformed = transformed - transformed.mean(axis=1, keepdims=True)
    pca = PCA(n_components=2, svd_solver="full")
    scores = pca.fit_transform(transformed)
    loadings = pca.components_.T
    _, loadings = orient_old_positive(scores, loadings, ages)
    return loadings[:, 0]


def main() -> None:
    cache_dir = Path(tempfile.gettempdir()) / "pflog_matplotlib_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    adata = ad.read_h5ad(DEFAULT_INPUT)
    counts = as_csr(adata.X)
    depths = np.asarray(counts.sum(axis=1)).ravel()
    alpha, sbar = estimate_alpha(counts)
    y0_opt = 1.0 / (4.0 * alpha)
    reference = {
        "PFlog optimum": y0_opt,
        "CPM": sbar / 1_000_000.0,
        "log1pPF": 1.0,
        "CP10k": sbar / 10_000.0,
    }

    pf_counts = counts.multiply(sbar / depths[:, None]).toarray()
    ages = adata.obs["age_label"].astype(str).to_numpy()
    de = pd.read_csv(DEFAULT_DE).dropna(subset=["logFC"]).copy()
    gene_index = pd.Series(np.arange(adata.n_vars), index=adata.var_names.astype(str))
    idx = gene_index.loc[de["gene"].astype(str)].to_numpy()
    logfc = de["logFC"].to_numpy(dtype=float)

    lo = min(reference.values()) / 8.0
    hi = max(reference.values()) * 8.0
    pseudocounts = np.unique(np.r_[np.logspace(np.log10(lo), np.log10(hi), 90), list(reference.values())])

    rows = []
    for y0 in pseudocounts:
        for label, center in [("logPF", False), ("PFlog", True)]:
            loading = pc1_loadings_for_pseudocount(pf_counts, ages, float(y0), clr_center=center)
            r = pearsonr(loading[idx], logfc).statistic
            rows.append({"method": label, "pseudocount": float(y0), "r2": float(r * r)})
    sweep = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    table_out = DEFAULT_OUT.with_suffix(".tsv")
    sweep.to_csv(table_out, sep="\t", index=False)

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
    colors = {"logPF": "#4d4d4d", "PFlog": "#e41a1c"}
    for method, group in sweep.groupby("method", sort=False):
        group = group.sort_values("pseudocount")
        ax.plot(group["pseudocount"], group["r2"], color=colors[method], lw=1.8, label=method)

    marker_styles = {
        "PFlog optimum": dict(color="#e41a1c", linestyle="-", lw=1.0),
        "CPM": dict(color="#2c7fb8", linestyle="--", lw=0.9),
        "log1pPF": dict(color="#238b45", linestyle="--", lw=0.9),
        "CP10k": dict(color="#756bb1", linestyle="--", lw=0.9),
    }
    label_y = {"PFlog optimum": 0.05, "CPM": 0.14, "log1pPF": 0.36, "CP10k": 0.32}
    for name, value in reference.items():
        style = marker_styles[name]
        ax.axvline(value, **style)
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
    ax.set_ylim(0, min(1.0, max(0.92, sweep["r2"].max() * 1.08)))
    ax.set_xlabel("count-scale pseudocount")
    ax.set_ylabel(r"$R^2$(PC1 loading, edgePython logFC)")
    ax.legend(frameon=False, loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.7, length=3)

    fig.savefig(DEFAULT_OUT, bbox_inches="tight")
    fig.savefig(DEFAULT_OUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    summary = {
        "alpha": alpha,
        "mean_depth": sbar,
        "reference_pseudocounts": reference,
        "output_pdf": str(DEFAULT_OUT),
        "output_png": str(DEFAULT_OUT.with_suffix(".png")),
        "output_table": str(table_out),
    }
    (DEFAULT_OUT.with_suffix(".json")).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
