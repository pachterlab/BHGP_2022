#!/usr/bin/env python3
"""Across-dataset 3-panel summary boxplot of per-method metrics.

Walks every <data_root>/<DS>/subset_genes/<DS>_subset_genes_metrics.json,
filters datasets by minimum avg_per_cell (using a reference dataset as the
threshold), and produces a 3-panel boxplot over methods of:
    A) Coefficient of variation of gene variance  (variance stabilization)
    B) Pearson r^2 of cell depth vs raw            (depth normalization)
    C) 1 - |mean per-cell Spearman r|              (monotonicity, lower=better)

The reference dataset's points are overlaid as red dots.

Usage:
    python plot_summary.py <data_root> <out_prefix> [reference_ds]

Defaults: reference_ds = angelidis_2019.
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 15
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LABELS = [
    "raw",
    "PF",
    "sqrt",
    "log1p",
    "log1pCP10k",
    "log1pCPM",
    "scalelog1pCP10k",
    "sctransform",
    "log1pPF",
    "PFlogPF (shift. CLR)",
]
TOP_KEYS = [
    "ncells", "ngenes", "density", "total_counts",
    "avg_per_cell", "avg_per_gene", "min_cell", "max_cell",
    "nvals", "overdispersion",
]

cividis = matplotlib.colormaps["cividis"]
COLORS = {"cell": cividis(0.01), "gene": cividis(0.5), "mono": cividis(0.99)}


def load_metrics(data_root):
    rows = []
    pattern = os.path.join(data_root, "*", "subset_genes", "*_subset_genes_metrics.json")
    for fn in sorted(glob.glob(pattern)):
        ds = os.path.basename(os.path.dirname(os.path.dirname(fn)))
        with open(fn) as f:
            m = json.load(f)
        top = {k: m.get(k) for k in TOP_KEYS}
        for method in LABELS:
            mm = m.get(method)
            if not isinstance(mm, dict):
                continue
            rows.append({
                "ds": ds,
                "method": method,
                **top,
                "cov_gene": mm.get("cov_gene"),
                "cov_cell": mm.get("cov_cell"),
                "r2_depth": mm.get("r2_depth"),
                "r_mono": mm.get("r_mono"),
            })
    return pd.DataFrame(rows)


def main(data_root, out_prefix, reference_ds="angelidis_2019"):
    df = load_metrics(data_root)
    # Drop datasets that have any NaN metric in any of the LABELS we will plot.
    bad = df.groupby("ds")[["cov_gene", "r2_depth", "r_mono"]].apply(lambda x: x.isna().any().any())
    keep_ds_all = bad[~bad].index
    df = df[df["ds"].isin(keep_ds_all)]
    n_total = df["ds"].nunique()

    # Min-depth filter at the reference dataset's avg_per_cell.
    ref_rows = df[(df["method"] == "raw") & (df["ds"] == reference_ds)]
    if ref_rows.empty:
        print(f"WARNING: reference dataset '{reference_ds}' missing — no depth filter applied", file=sys.stderr)
        df_pass = df
    else:
        min_apc = float(ref_rows["avg_per_cell"].iloc[0])
        pass_ds = df[(df["method"] == "raw") & (df["avg_per_cell"] >= min_apc)]["ds"].unique()
        df_pass = df[df["ds"].isin(pass_ds)]
    n_pass = df_pass["ds"].nunique()
    print(f"{n_total} datasets after NaN drop, {n_pass} pass the avg_per_cell >= {reference_ds} filter")

    # mono metric: lower = more monotonic
    df_pass = df_pass.copy()
    df_pass["mono_metric"] = 1 - df_pass["r_mono"].abs()

    fig, axs = plt.subplots(figsize=(7.5, 12), nrows=3)
    fig.subplots_adjust(hspace=0)
    fig.suptitle(
        f"{n_pass:,} of {n_total:,} pass-filter datasets (subset_genes)",
        y=0.905,
    )

    panels = [
        ("cov_gene",    "gene", "Coefficient of variation\ngene variance"),
        ("r2_depth",    "cell", "Pearson $r^2$ cell depth"),
        ("mono_metric", "mono", r"$1 - |$mean Spearman $r|$"),
    ]

    gb = df_pass.groupby("method")
    unique_size = 9 ** 2

    ref_in_pass = reference_ds in df_pass["ds"].unique()
    for ax, (met, color_key, ylabel) in zip(axs, panels):
        per_method = gb[met].apply(lambda s: s.values).reindex(LABELS)
        # Drop methods that may be missing entirely.
        present_labels = [m for m in LABELS if isinstance(per_method.loc[m], np.ndarray)]
        data = [per_method.loc[m] for m in present_labels]
        positions = np.arange(len(present_labels))

        bp = ax.boxplot(
            data, positions=positions, patch_artist=True, notch=True,
            boxprops=dict(facecolor=COLORS[color_key]),
            medianprops=dict(color="k"),
            flierprops=dict(markerfacecolor=COLORS[color_key], markeredgecolor="k", markersize=4),
            widths=0.6,
        )

        means = np.array([np.nanmean(d) for d in data])
        ax.scatter(positions, means, color="k", s=unique_size * 17, linewidth=3, marker="_", zorder=10)
        mean_dot = ax.scatter(positions, means, facecolor="white", edgecolor="k", s=unique_size, zorder=10, label="mean")

        if ref_in_pass:
            ref_y = []
            for m in present_labels:
                vv = df_pass[(df_pass["method"] == m) & (df_pass["ds"] == reference_ds)][met].values
                ref_y.append(vv[0] if len(vv) else np.nan)
            ref_y = np.array(ref_y)
            ref_dot = ax.scatter(
                positions, ref_y, facecolor="#DC3220", edgecolor="k", s=unique_size,
                zorder=11, label=reference_ds,
            )

        ax.set_ylabel(ylabel, labelpad=12)
        ax.yaxis.set_label_coords(-0.10, 0.5)
        if met == "cov_gene":
            ax.set(yscale="symlog", ylim=(-1, ax.get_ylim()[1] if ax.get_ylim()[1] > 1 else 250))
        elif met == "r2_depth":
            ax.set_ylim(-0.05, 1.05)
        else:
            ax.set_ylim(-0.05, 1.05)

        if ax is not axs[-1]:
            ax.set_xticks(positions)
            ax.set_xticklabels([])
        else:
            ax.set_xticks(positions)
            ax.set_xticklabels(present_labels, rotation=45, ha="right")

        ax.legend(loc="best", prop={"size": 12})

    out_dir = os.path.dirname(out_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(f"{out_prefix}.pdf", facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(f"{out_prefix}.png", facecolor="white", bbox_inches="tight", dpi=200)
    print(f"wrote {out_prefix}.pdf + .png")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    data_root = sys.argv[1]
    out_prefix = sys.argv[2]
    ref = sys.argv[3] if len(sys.argv) > 3 else "angelidis_2019"
    main(data_root, out_prefix, ref)
