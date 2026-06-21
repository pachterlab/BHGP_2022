#!/usr/bin/env python3
"""Generate manuscript Figure 2 from Angelidis lung pseudobulk counts.

The input AnnData contains raw UMI counts summed by mouse sample. The script
compares standard log1p(CP10K) normalization with PFlog normalization, then
compares PC1 gene loadings with edgePython old-vs-young differential expression.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

import scclr


matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
AEH_BENCHMARK_DIR = HERE.parents[1]
FIGURE_OUT = AEH_BENCHMARK_DIR / "output"

DEFAULT_INPUT = HERE / "angelidis_lung_pseudobulk.h5ad"
DEFAULT_FIGURE = FIGURE_OUT / "angelidis_pca_and_pc1_loadings_six_panel.pdf"
DEFAULT_PSEUDOCOUNT_SWEEP = FIGURE_OUT / "angelidis_pseudobulk_pseudocount_r2_sweep.tsv"
DEFAULT_PSEUDOCOUNT_SUMMARY = FIGURE_OUT / "angelidis_pseudobulk_pseudocount_r2_sweep.json"

AGE_COLORS = {"Young": "#2c7fb8", "Old": "#d95f5f"}
AGE_MARKERS = {"Young": "o", "Old": "^"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--pseudocount-sweep", type=Path, default=DEFAULT_PSEUDOCOUNT_SWEEP)
    parser.add_argument("--pseudocount-summary", type=Path, default=DEFAULT_PSEUDOCOUNT_SUMMARY)
    parser.add_argument("--n-comps", type=int, default=5)
    parser.add_argument(
        "--comparison",
        choices=("cp10k", "same-k"),
        default="cp10k",
        help="Comparison transform for panels a/c. 'same-k' uses log(1 + K x/s) with the PFlog-estimated K.",
    )
    parser.add_argument(
        "--edgepython-path",
        type=Path,
        default=os.environ.get("EDGEPYTHON_PATH"),
        help="Optional path containing the edgepython package.",
    )
    return parser.parse_args()


def as_csr(x) -> sp.csr_matrix:
    return x.tocsr() if sp.issparse(x) else sp.csr_matrix(np.asarray(x, dtype=np.float64))


def orient_old_positive(scores: np.ndarray, loadings: np.ndarray, obs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    old = obs["age_label"].astype(str).to_numpy() == "Old"
    young = obs["age_label"].astype(str).to_numpy() == "Young"
    if scores[old, 0].mean() < scores[young, 0].mean():
        scores = scores.copy()
        loadings = loadings.copy()
        scores[:, 0] *= -1.0
        loadings[:, 0] *= -1.0
    return scores, loadings


def run_pflog_pca(adata: ad.AnnData, n_comps: int) -> ad.AnnData:
    work = adata.copy()
    scclr.pp.pflog(work, target="auto")
    scclr.tl.pca(work, n_comps=n_comps, ncv=min(12, work.n_obs - 1), seed=0)
    scores, loadings = orient_old_positive(work.obsm["X_pca"], work.varm["PCs"], work.obs)
    work.obsm["X_pca"] = scores
    work.varm["PCs"] = loadings
    return work


def run_logcp10k_pca(adata: ad.AnnData, n_comps: int) -> ad.AnnData:
    work = adata.copy()
    counts = as_csr(work.X)
    depths = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(depths <= 0):
        raise ValueError("all samples must have positive total counts")

    normalized = sp.diags(10000.0 / depths) @ counts
    normalized.data = np.log1p(normalized.data)
    work.layers["logcp10k"] = normalized.tocsr()

    pca = PCA(n_components=n_comps, svd_solver="full")
    scores = pca.fit_transform(normalized.toarray())
    loadings = pca.components_.T
    scores, loadings = orient_old_positive(scores, loadings, work.obs)
    work.obsm["X_pca"] = np.ascontiguousarray(scores)
    work.varm["PCs"] = np.ascontiguousarray(loadings)
    work.uns["pca"] = {
        "variance": pca.explained_variance_,
        "variance_ratio": pca.explained_variance_ratio_,
    }
    work.uns["method_label"] = "log1p(CP10K)"
    work.uns["loading_label"] = "log1p(CP10K) PC1 loading"
    return work


def run_log_same_k_pca(adata: ad.AnnData, n_comps: int, k_scale: float) -> ad.AnnData:
    work = adata.copy()
    counts = as_csr(work.X)
    depths = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(depths <= 0):
        raise ValueError("all samples must have positive total counts")

    normalized = sp.diags(float(k_scale) / depths) @ counts
    normalized.data = np.log1p(normalized.data)
    work.layers["log_same_k"] = normalized.tocsr()

    pca = PCA(n_components=n_comps, svd_solver="full")
    scores = pca.fit_transform(normalized.toarray())
    loadings = pca.components_.T
    scores, loadings = orient_old_positive(scores, loadings, work.obs)
    work.obsm["X_pca"] = np.ascontiguousarray(scores)
    work.varm["PCs"] = np.ascontiguousarray(loadings)
    work.uns["pca"] = {
        "variance": pca.explained_variance_,
        "variance_ratio": pca.explained_variance_ratio_,
    }
    work.uns["method_label"] = f"log(1 + K x/s), K={k_scale:.3g}"
    work.uns["loading_label"] = "logPF same-K PC1 loading"
    return work


def run_edgepython_de(adata: ad.AnnData, edgepython_path: Path | None) -> pd.DataFrame:
    if edgepython_path is not None:
        sys.path.insert(0, str(edgepython_path))

    try:
        import edgepython as ep
    except Exception as exc:
        cached = HERE / "pc1_loadings_with_edgepython_de.csv"
        if cached.exists():
            cached_de = pd.read_csv(cached)
            loading_cols = {"pflog_pc1_loading", "logcp10k_pc1_loading"}
            de_cols = [col for col in cached_de.columns if col not in loading_cols]
            print(f"edgePython unavailable ({exc}); reusing DE columns from {cached}", file=sys.stderr)
            return cached_de[de_cols].copy()
        raise

    counts = as_csr(adata.X).toarray().T.astype(np.float64, copy=False)
    genes = adata.var_names.astype(str).to_numpy()
    condition = (adata.obs["age_label"].astype(str).to_numpy() == "Old").astype(int)

    y = ep.make_dgelist(counts, group=condition)
    keep = ep.filter_by_expr(y, group=condition)
    y = ep.make_dgelist(counts[keep], group=condition)
    y = ep.calc_norm_factors(y)
    design = np.column_stack([np.ones(len(condition)), condition])

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, module=r"edgepython\..*")
        y = ep.estimate_disp(y, design)
        fit = ep.glm_ql_fit(y, design)
        qlf = ep.glm_ql_ftest(fit, coef=1)

    table = ep.top_tags(qlf, n=int(keep.sum()))["table"].copy()
    filtered_genes = genes[keep]
    table.insert(0, "gene", filtered_genes[table.index.to_numpy()])
    return table.reset_index(drop=True)


def write_pca_coordinates(adata_log: ad.AnnData, adata_pf: ad.AnnData) -> pd.DataFrame:
    rows = []
    for method, work in [(adata_log.uns.get("method_label", "log1p(CP10K)"), adata_log), ("PFlog", adata_pf)]:
        for i, sample in enumerate(work.obs_names):
            rows.append(
                {
                    "method": method,
                    "sample": sample,
                    "PC1": work.obsm["X_pca"][i, 0],
                    "PC2": work.obsm["X_pca"][i, 1],
                    "age": work.obs["age"].iloc[i],
                    "age_label": work.obs["age_label"].iloc[i],
                    "n_cells": int(work.obs["n_cells"].iloc[i]),
                    "total_umi": float(work.obs["total_umi"].iloc[i]),
                }
            )
    coords = pd.DataFrame(rows)
    coords.to_csv(HERE / "angelidis_pseudobulk_pca_coordinates.csv", index=False)
    return coords


def merge_loadings_with_de(adata: ad.AnnData, adata_log: ad.AnnData, adata_pf: ad.AnnData, de: pd.DataFrame) -> pd.DataFrame:
    loadings = pd.DataFrame(
        {
            "gene": adata.var_names.astype(str),
            "pflog_pc1_loading": adata_pf.varm["PCs"][:, 0],
            "logcp10k_pc1_loading": adata_log.varm["PCs"][:, 0],
        }
    )
    merged = loadings.merge(de, on="gene", how="left")
    merged.to_csv(HERE / "pc1_loadings_with_edgepython_de.csv", index=False)
    return merged


def panel_label(ax, label: str) -> None:
    ax.text(-0.14, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top", ha="left")


def style_legend(legend) -> None:
    frame = legend.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("#4d4d4d")
    frame.set_linewidth(0.6)
    frame.set_alpha(0.92)


def plot_pca_panel(
    ax,
    work: ad.AnnData,
    panel: str,
    method_label: str,
    legend_loc: str = "lower right",
    legend_bbox: tuple[float, float] | None = None,
    borderaxespad: float = 0.5,
) -> None:
    scores = work.obsm["X_pca"]
    var = work.uns["pca"]["variance_ratio"]
    for label in ["Young", "Old"]:
        mask = work.obs["age_label"].astype(str).to_numpy() == label
        ax.scatter(
            scores[mask, 0],
            scores[mask, 1],
            s=40,
            marker=AGE_MARKERS[label],
            c=AGE_COLORS[label],
            edgecolor="black",
            linewidth=0.45,
            alpha=0.95,
            label=label,
        )
        for sample, x, y in zip(work.obs_names[mask], scores[mask, 0], scores[mask, 1]):
            ax.annotate(sample, (x, y), xytext=(3, 2), textcoords="offset points", fontsize=5.8)

    ax.axhline(0, color="#d9d9d9", linewidth=0.7, zorder=0)
    ax.axvline(0, color="#d9d9d9", linewidth=0.7, zorder=0)
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}% variance)")
    ax.text(0.05, 0.94, method_label, transform=ax.transAxes, fontsize=8, va="top", ha="left")
    legend_kwargs = dict(
        frameon=True,
        fancybox=False,
        loc=legend_loc,
        handletextpad=0.4,
        borderpad=0.25,
        borderaxespad=borderaxespad,
    )
    if legend_bbox is not None:
        legend_kwargs["bbox_to_anchor"] = legend_bbox
    legend = ax.legend(**legend_kwargs)
    style_legend(legend)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.7, length=3)
    panel_label(ax, panel)


def plot_depth_panel(ax, adata: ad.AnnData, panel: str) -> None:
    obs = adata.obs.copy().reset_index(drop=True)
    obs["depth_millions"] = obs["total_umi"].astype(float) / 1_000_000.0
    obs["age_label"] = obs["age_label"].astype(str)

    groups = ["Young", "Old"]
    data = [obs.loc[obs["age_label"] == label, "depth_millions"].to_numpy() for label in groups]
    positions = np.arange(len(groups))
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 0.9},
        boxprops={"linewidth": 0.75},
        whiskerprops={"linewidth": 0.75},
        capprops={"linewidth": 0.75},
    )
    for patch, label in zip(bp["boxes"], groups):
        patch.set_facecolor(AGE_COLORS[label])
        patch.set_alpha(0.22)
        patch.set_edgecolor("black")

    rng = np.random.default_rng(0)
    for i, label in enumerate(groups):
        vals = obs.loc[obs["age_label"] == label, "depth_millions"].to_numpy()
        jitter = rng.uniform(-0.08, 0.08, size=vals.shape[0])
        ax.scatter(
            np.full(vals.shape[0], positions[i]) + jitter,
            vals,
            s=28,
            color=AGE_COLORS[label],
            edgecolor="black",
            linewidth=0.45,
            alpha=0.95,
            zorder=3,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(groups)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("total UMI (millions)")
    ax.set_xlabel("age group")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.7, length=3)
    panel_label(ax, panel)


def plot_pseudocount_panel(ax, sweep_path: Path, summary_path: Path, panel: str) -> None:
    if not sweep_path.exists():
        raise FileNotFoundError(f"pseudocount sweep table not found: {sweep_path}")
    if not summary_path.exists():
        raise FileNotFoundError(f"pseudocount sweep summary not found: {summary_path}")

    sweep = pd.read_csv(sweep_path, sep="\t")
    summary = json.loads(summary_path.read_text())
    reference = summary["reference_pseudocounts"]

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
    label_y = {"PFlog optimum": 0.48, "CPM": 0.05, "log1pPF": 0.20, "CP10k": 0.32}
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
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.7, length=3)
    panel_label(ax, panel)


def plot_figure(
    adata_log: ad.AnnData,
    adata_pf: ad.AnnData,
    merged: pd.DataFrame,
    output: Path,
    pseudocount_sweep: Path,
    pseudocount_summary: Path,
) -> dict:
    cache_dir = Path(tempfile.gettempdir()) / "pflog_matplotlib_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    tested = merged.dropna(subset=["logFC", "FDR"]).copy()
    sig = tested["FDR"] < 0.05

    matplotlib.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
        }
    )

    fig, axes = plt.subplots(3, 2, figsize=(7.1, 8.2), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.03, h_pad=0.03, wspace=0.08, hspace=0.08)

    comparison_label = adata_log.uns.get("method_label", "log1p(CP10K)")
    comparison_loading_label = adata_log.uns.get("loading_label", "log1p(CP10K) PC1 loading")

    plot_pca_panel(axes[0, 0], adata_log, "a", f"{comparison_label} PCA", legend_loc="upper right", borderaxespad=0.18)
    plot_pca_panel(axes[0, 1], adata_pf, "b", "PFlog PCA", legend_loc="lower right", legend_bbox=(0.83, 0.03))

    scatter_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#2c7fb8",
            markeredgecolor="none",
            markersize=4.5,
            alpha=0.75,
            label="FDR < 0.05",
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#bdbdbd",
            markeredgecolor="none",
            markersize=4.5,
            alpha=0.55,
            label="not significant",
        ),
    ]

    correlations = {}
    for ax, (panel, col, xlabel) in zip(
        axes[1, :],
        [
            ("c", "logcp10k_pc1_loading", comparison_loading_label),
            ("d", "pflog_pc1_loading", "PFlog PC1 loading"),
        ],
    ):
        r = pearsonr(tested[col], tested["logFC"]).statistic
        correlations[col] = float(r)
        ax.scatter(
            tested.loc[~sig, col],
            tested.loc[~sig, "logFC"],
            s=5.5,
            color="#bdbdbd",
            alpha=0.42,
            linewidths=0,
            rasterized=True,
        )
        ax.scatter(
            tested.loc[sig, col],
            tested.loc[sig, "logFC"],
            s=6.5,
            color="#2c7fb8",
            alpha=0.62,
            linewidths=0,
            rasterized=True,
        )
        ax.axhline(0, color="#d9d9d9", linewidth=0.7, zorder=0)
        ax.axvline(0, color="#d9d9d9", linewidth=0.7, zorder=0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("edgePython logFC (old vs young)")
        ax.text(0.05, 0.93, f"$R^2$ = {r ** 2:.2f}", transform=ax.transAxes, fontsize=8, va="top", ha="left")
        ax.legend(handles=scatter_handles, frameon=False, loc="lower right", handletextpad=0.4, borderpad=0.2)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=0.7, length=3)
        panel_label(ax, panel)

    plot_depth_panel(axes[2, 0], adata_pf, "e")
    plot_pseudocount_panel(axes[2, 1], pseudocount_sweep, pseudocount_summary, "f")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")

    return {
        "n_genes_tested": int(len(tested)),
        "n_fdr_0_05": int(sig.sum()),
        "pflog_k": float(adata_pf.uns["pflog"]["k"]),
        "pflog_pc1_variance_ratio": float(adata_pf.uns["pca"]["variance_ratio"][0]),
        "pflog_pc2_variance_ratio": float(adata_pf.uns["pca"]["variance_ratio"][1]),
        "logcp10k_pc1_variance_ratio": float(adata_log.uns["pca"]["variance_ratio"][0]),
        "logcp10k_pc2_variance_ratio": float(adata_log.uns["pca"]["variance_ratio"][1]),
        "pearson_pflog_loading_vs_logFC": correlations["pflog_pc1_loading"],
        "pearson_logcp10k_loading_vs_logFC": correlations["logcp10k_pc1_loading"],
    }


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    adata_pf = run_pflog_pca(adata, args.n_comps)
    if args.comparison == "same-k":
        adata_log = run_log_same_k_pca(adata, args.n_comps, float(adata_pf.uns["pflog"]["k"]))
    else:
        adata_log = run_logcp10k_pca(adata, args.n_comps)
    write_pca_coordinates(adata_log, adata_pf)
    de = run_edgepython_de(adata, args.edgepython_path)
    merged = merge_loadings_with_de(adata, adata_log, adata_pf, de)
    summary = plot_figure(adata_log, adata_pf, merged, args.output, args.pseudocount_sweep, args.pseudocount_summary)
    summary.update(
        {
            "n_samples": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_cells": int(adata.obs["n_cells"].sum()),
            "n_young": int((adata.obs["age_label"].astype(str) == "Young").sum()),
            "n_old": int((adata.obs["age_label"].astype(str) == "Old").sum()),
        }
    )
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"wrote: {args.output}")


if __name__ == "__main__":
    main()
