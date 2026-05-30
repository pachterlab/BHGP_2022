#!/usr/bin/env python3
"""3-panel bar plot of cell-type metrics produced by metrics_celltype.py.

Panels:
  (a) PC1 loadings entropy (fraction of max)
  (b) # false-positive DE genes (top-500 vs bottom-500 by raw depth, q<0.01)
  (c) Mean | pairwise Spearman r | within cell type

Bars are colored cividis by default. PFlogPF (shift. CLR) is colored red so
it matches the AE&H benchmark figure's compositional family.

Usage:
    python plot_celltype_bar.py <metrics.json> <out_prefix>
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
    "raw", "PF", "sqrt", "log1p", "log1pCP10k", "log1pCPM",
    "scalelog1pCP10k", "sctransform", "log1pPF",
    "PFlogPF (shift. CLR)",
]
CLR_METHOD = "PFlogPF (shift. CLR)"
CLR_COLOR  = "#E41A1C"   # matches the Comp. family color in main_benchmark_fig

cividis = matplotlib.colormaps["cividis"]
PANEL_COLORS = {
    "pc1":      cividis(0.5),    # gene-axis tone
    "fp":       cividis(0.01),   # depth/cell tone
    "spearman": cividis(0.99),   # monotonicity/yellow tone
}


def main(metrics_json, out_prefix):
    with open(metrics_json) as f:
        d = json.load(f)
    ds, ct = d.get("dataset", "?"), d.get("celltype", "?")
    n_cells, n_genes = d.get("n_cells", 0), d.get("n_genes", 0)

    methods_present = [m for m in METHODS if m in d["methods"]]
    # Panel (a) is plotted as 1 - fraction so depth-confounded (loadings
    # concentrated in few genes) reads as a tall bar, matching the original
    # manuscript figure.
    pc1 = np.array([1 - d["methods"][m]["pc1_entropy_frac"] for m in methods_present])
    fp  = np.array([d["methods"][m]["fp_de_genes"]     for m in methods_present])
    # Panel (c) is plotted as 1 - |mean Spearman| (low = cells correlate well,
    # high = transform crushed cell-to-cell variability).
    sp  = np.array([1 - abs(d["methods"][m]["mean_abs_spearman"]) for m in methods_present])

    x = np.arange(len(methods_present))

    fig, axs = plt.subplots(nrows=3, figsize=(12, 11))
    fig.subplots_adjust(hspace=0.05)
    title = ct.replace("_", " ")
    fig.suptitle(f"{title}  ({n_cells:,} cells × {n_genes:,} genes, {ds})",
                 y=0.905)

    panels = [
        ("pc1",      pc1, "1 - Fraction of max entropy on PC1 loadings"),
        ("fp",       fp,  "Number of DE genes\n(highest 500 vs lowest 500 raw depth)"),
        ("spearman", sp,  r"$1-|$mean within-celltype pairwise Spearman $r|$"),
    ]

    for ax, (key, y, ylabel) in zip(axs, panels):
        face = [CLR_COLOR if m == CLR_METHOD else PANEL_COLORS[key]
                for m in methods_present]
        ax.bar(x, y, width=0.75, facecolor=face, edgecolor="k")
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        if ax is axs[-1]:
            ax.set_xticklabels(methods_present, rotation=45, ha="right")
        else:
            ax.set_xticklabels([])
        # Show numerical value on top of each bar (helps when ranges differ).
        for xi, yi in zip(x, y):
            ax.text(xi, yi, f"{yi:,.2g}" if key != "fp" else f"{int(yi)}",
                    ha="center", va="bottom", fontsize=9)

    # Tight y limits per panel
    axs[0].set_ylim(0, max(1.05, np.nanmax(pc1) * 1.10))
    # Linear scale — symlog squashes the differences between methods that the
    # bar plot is meant to highlight (e.g. PF/PFlogPF (shift. CLR) being much
    # lower than the rest).
    if fp.max() > 0:
        axs[1].set_ylim(0, fp.max() * 1.10)
    axs[2].set_ylim(0, 1.05)

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    fig.savefig(f"{out_prefix}.pdf", facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(f"{out_prefix}.png", facecolor="white", bbox_inches="tight", dpi=200)
    print(f"wrote {out_prefix}.pdf + .png")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr); sys.exit(1)
    main(sys.argv[1], sys.argv[2])
