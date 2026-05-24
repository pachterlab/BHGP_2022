#!/usr/bin/env python3
"""3-panel bar plot of per-method metrics for a single dataset.

Reads the per-dataset metrics JSON produced by metrics_methods_*.py and
draws three bar panels:
    A) Coefficient of variation of gene variances  (variance stabilization)
    B) Pearson r2 of cell depth vs raw             (depth normalization)
    C) Mean per-cell Spearman r vs raw             (monotonicity)

Usage:
    python plot_bar.py <metrics.json> <out_prefix>
"""
import json
import os
import sys

import matplotlib
matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 14
import matplotlib.pyplot as plt
import numpy as np

METHODS = [
    "raw",
    "PF",
    "sqrt",
    "log1p",
    "log1pCP10k",
    "log1pCPM",
    "scalelog1pCP10k",
    "sctransform",
    "log1pPF",
    "PFlog1pPF (CLR)",
]

cividis = matplotlib.colormaps["cividis"]
COLORS = {"gene": cividis(0.5), "cell": cividis(0.01), "mono": cividis(0.99)}


def main(metrics_json, out_prefix):
    with open(metrics_json) as f:
        metrics = json.load(f)

    ds = os.path.basename(metrics_json).split("_subset_genes")[0].split("_all_genes")[0]
    ncells = metrics.get("ncells", "?")
    ngenes = metrics.get("ngenes", "?")

    def get(method, key):
        v = metrics.get(method, {}).get(key)
        return float("nan") if v is None else float(v)

    cov_gene = np.array([get(m, "cov_gene") for m in METHODS])
    r2_depth = np.array([get(m, "r2_depth") for m in METHODS])
    r_mono = np.array([get(m, "r_mono") for m in METHODS])
    mono_metric = 1 - np.abs(r_mono)

    x_idx = np.arange(len(METHODS))
    fig, axs = plt.subplots(nrows=3, figsize=(10, 4.5 * 3))
    title = f"{ds} ({ncells:,} cells x {ngenes:,} genes)" if isinstance(ncells, int) else str(ds)
    fig.suptitle(title, y=0.92)
    fig.subplots_adjust(hspace=0.55)

    ax = axs[0]
    ax.bar(x_idx, cov_gene, width=0.75, facecolor=COLORS["gene"], edgecolor="k")
    ax.set(yscale="symlog", ylabel="Coefficient of variation\ngene variance", xticks=[])
    ax.set_title("(A) Variance stabilization", loc="left")

    ax = axs[1]
    ax.bar(x_idx, r2_depth, width=0.75, facecolor=COLORS["cell"], edgecolor="k")
    ax.set(ylabel="Pearson $r^2$ cell depth", ylim=(0, 1.05), xticks=[])
    ax.set_title("(B) Depth normalization", loc="left")

    ax = axs[2]
    ax.bar(x_idx, mono_metric, width=0.75, facecolor=COLORS["mono"], edgecolor="k")
    ax.set(ylabel=r"$1 - |$mean Spearman $r|$", ylim=(0, 1.05), xticks=x_idx)
    ax.set_xticklabels(METHODS, rotation=45, ha="right")
    ax.set_title("(C) Monotonicity (lower = more monotonic)", loc="left")

    fig.savefig(f"{out_prefix}.pdf", facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(f"{out_prefix}.png", facecolor="white", bbox_inches="tight", dpi=200)
    print(f"wrote {out_prefix}.pdf + .png")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
