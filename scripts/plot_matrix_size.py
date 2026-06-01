#!/usr/bin/env python3
"""Panel (b) of the speed/memory supplement: matrix size vs number of cells.

Reproduces the scatter of in-memory matrix size [GB] against the number of
cells for the sparse raw matrix (green) and the dense sctransform matrix
(red), over every dataset in the benchmark. Companion to panel (a), which is
a reproduction of Fig 3c from Choudhary & Satija 2022.

Reads the tab-separated size tables (ds, normalization, ncells, nbytes)
produced during the benchmark:
    <data_root>/size_raw.txt
    <data_root>/size_sctransform.txt

Usage:
    python plot_matrix_size.py <data_root> <out_prefix>
"""
import os
import sys

import matplotlib
matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 15
import matplotlib.pyplot as plt
import pandas as pd

COLS = ["ds", "normalization", "ncells", "nbytes"]


def load(data_root, norm):
    fn = os.path.join(data_root, f"size_{norm}.txt")
    return pd.read_csv(fn, header=None, sep="\t", names=COLS)


def main(data_root, out_prefix):
    raw = load(data_root, "raw")
    sct = load(data_root, "sctransform")

    fig, ax = plt.subplots(figsize=(5, 5))

    y_raw = raw["nbytes"] / 1e9
    ax.scatter(raw["ncells"], y_raw, edgecolor="k", facecolor="green",
               label=f"raw ({y_raw.max():,.2f} GB max)")

    y_sct = sct["nbytes"] / 1e9
    ax.scatter(sct["ncells"], y_sct, edgecolor="k", facecolor="red",
               label="sctransform")

    ax.set(xlabel="Number of cells", ylabel="Matrix size [GB]", xscale="log")
    ax.legend()

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    fig.savefig(f"{out_prefix}.pdf", facecolor="white", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}.png", facecolor="white", dpi=300, bbox_inches="tight")
    print(f"wrote {out_prefix}.pdf + .png "
          f"({len(raw)} datasets, raw max {y_raw.max():.2f} GB, "
          f"sctransform max {y_sct.max():.2f} GB)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
