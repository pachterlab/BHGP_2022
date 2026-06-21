#!/usr/bin/env python3
"""3-panel bar plot for Fig. 1a.

Panels:
  (a) PC1 loadings entropy (fraction of max)
  (b) mean significant abundance-matched log-ratio balance pairs
  (c) Mean | pairwise Spearman r | within cell type

Bars are colored cividis by default. PFlog (shift. CLR) is colored red so
it matches the AE&H benchmark figure's compositional family.

Usage:
    python plot_celltype_bar.py <metrics.json> <out_prefix> \
        --balance-summary analysis/figures/angelidis_type2_3m_many_balance_de_summary.tsv
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 14
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

METHODS = [
    "raw", "PF", "sqrt", "log1p", "log1pCP10k", "log1pCPM",
    "sctransform", "log1pPF",
    "PFlog (shift. CLR)",
]
CLR_METHOD = "PFlog (shift. CLR)"
CLR_COLOR  = "#E41A1C"   # matches the Comp. family color in main_benchmark_fig
DISPLAY = {"sctransform": "sctransform v2"}   # shown text only; data key unchanged
BALANCE_KEY = {"PFlog (shift. CLR)": "PFlog"}

cividis = matplotlib.colormaps["cividis"]
PANEL_COLORS = {
    "pc1":      cividis(0.5),    # gene-axis tone
    "fp":       cividis(0.01),   # depth/cell tone
    "spearman": cividis(0.99),   # monotonicity/yellow tone
}


def main(metrics_json, out_prefix, balance_summary):
    with open(metrics_json) as f:
        d = json.load(f)
    balance = pd.read_csv(balance_summary, sep="\t").set_index("method")

    methods_present = [m for m in METHODS if m in d["methods"]]
    # Panel (a) is plotted as 1 - fraction so depth-confounded (loadings
    # concentrated in few genes) reads as a tall bar, matching the original
    # manuscript figure.
    pc1 = np.array([1 - d["methods"][m]["pc1_entropy_frac"] for m in methods_present])
    balance_names = [BALANCE_KEY.get(m, m) for m in methods_present]
    de_pairs = np.array([balance.loc[m, "mean_significant_balances"] for m in balance_names])
    # Panel (c) is plotted as 1 - |mean Spearman| (low = cells correlate well,
    # high = transform crushed cell-to-cell variability).
    sp  = np.array([1 - abs(d["methods"][m]["mean_abs_spearman"]) for m in methods_present])

    x = np.arange(len(methods_present))

    fig, axs = plt.subplots(nrows=3, figsize=(7.5, 12))
    # Fixed absolute margins so the final canvas is exactly 7.5x12 in and
    # matches plot_summary.py's geometry pixel-for-pixel (needed for
    # side-by-side LaTeX placement).
    fig.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.18, hspace=0)

    panels = [
        ("pc1",      pc1, "1 - Fraction of max entropy on PC1 loadings"),
        ("fp",       de_pairs,  "Mean DE pairs"),
        ("spearman", sp,  r"1 - $|$ mean w/in-celltype pairwise Spearman $r$ $|$"),
    ]

    for ax, (key, y, ylabel) in zip(axs, panels):
        face = [CLR_COLOR if m == CLR_METHOD else PANEL_COLORS[key]
                for m in methods_present]
        ax.bar(x, y, width=0.75, color=face, edgecolor="k")
        ax.set_ylabel(ylabel, labelpad=12)
        # push the ylabel further left of the tick labels
        ax.yaxis.set_label_coords(-0.10, 0.5)
        ax.set_xticks(x)
        if ax is axs[-1]:
            ax.set_xticklabels([DISPLAY.get(m, m) for m in methods_present],
                               rotation=45, ha="right")
        else:
            ax.set_xticklabels([])
        # Show numerical value on top of each bar (helps when ranges differ).
        for xi, yi in zip(x, y):
            ax.text(xi, yi, f"{yi:,.2g}", ha="center", va="bottom", fontsize=9)

    # Tight y limits per panel
    axs[0].set_ylim(0, max(1.05, np.nanmax(pc1) * 1.10))
    # Linear scale — symlog squashes the differences between methods that the
    # bar plot is meant to highlight (e.g. PF/PFlog (shift. CLR) being much
    # lower than the rest).
    if de_pairs.max() > 0:
        axs[1].set_ylim(0, de_pairs.max() * 1.15)
    axs[2].set_ylim(0, 1.05)

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    # No bbox_inches="tight" — keep the full 7.5x12 in canvas so the PDF
    # matches summary_subset_genes.pdf in absolute dimensions.
    fig.savefig(f"{out_prefix}.pdf", facecolor="white", dpi=300)
    fig.savefig(f"{out_prefix}.png", facecolor="white", dpi=200)
    print(f"wrote {out_prefix}.pdf + .png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics_json")
    parser.add_argument("out_prefix")
    parser.add_argument(
        "--balance-summary",
        default="analysis/figures/angelidis_type2_3m_many_balance_de_summary.tsv",
    )
    args = parser.parse_args()
    main(args.metrics_json, args.out_prefix, args.balance_summary)
