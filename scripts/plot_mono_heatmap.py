#!/usr/bin/env python3
"""DE-gene-scramble / monotonicity heatmap (supplement Fig 8).

For one dataset (default angelidis_2019):
  1. Build a celltype x gene matrix per method by averaging the expression of
     every cell within a cell type.
  2. Pick the top-100 genes for each cell type by their average expression in
     PF (this defines a per-celltype gene set AND order, shared across panels).
  3. Rank those 100 genes within each cell type from lowest to highest with
     scipy.stats.rankdata(..., method="ordinal").
  4. Plot the matrix of ranks as a heatmap (cividis), one panel per method.

A truly monotonic transform preserves the PF gene ordering, so its heatmap is
a clean left->right gradient; rank-distorting transforms look scrambled.

Usage:
    python plot_mono_heatmap.py <dataset_dir> <out_prefix> [n_top]

<dataset_dir> must contain metadata_barcodes.txt.gz (last column 'celltype',
matrix-row order) and subset_genes/{pf,pf_log_pf}.mtx.gz +
subset_genes/{cp10k_log_scale,sctransform}.csv.gz.
"""
import gzip
import os
import sys

import matplotlib
matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 15
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread

# Display label -> (file under subset_genes/, is_dense_csv). Order = panel order
# (upper-left, upper-right, lower-left, lower-right).
PANELS = [
    ("raw",         "raw.mtx.gz",         False),
    ("log1pPF",     "pf_log.mtx.gz",      False),
    # PFlog = additive centered-log-ratio computed on the fly from raw via
    # norm_clr (PF to mean depth -> log1p -> per-cell centering); "__CLR__"
    # sentinel triggers that. NOT the multiplicative pf_log_pf.mtx.gz.
    ("PFlog",   "__CLR__",            True),
    ("sctransform", "sctransform.csv.gz", True),
]
REF = "raw"   # method whose per-celltype top-N order is used for every panel


def load_matrix(path, dense):
    if path.endswith(".bin"):
        n_cells, n_genes = map(int, open(path + ".shape").read().strip().split(","))
        return np.fromfile(path, dtype=np.float32).reshape(n_cells, n_genes)
    if dense:
        # pandas C parser is ~10x faster than np.loadtxt for these large csvs.
        return pd.read_csv(path, header=None, dtype=np.float32).values
    m = mmread(path)
    return np.asarray(m.todense(), dtype=np.float32)


def celltype_means(X, celltypes_sorted, assignments):
    """Return a (n_celltype x n_gene) matrix of per-celltype mean expression,
    rows ordered by celltypes_sorted."""
    df = pd.DataFrame(X, index=assignments)
    gg = df.groupby(df.index).mean()
    return gg.reindex(celltypes_sorted).values


def main(dataset_dir, out_prefix, n_top=100):
    meta = pd.read_csv(os.path.join(dataset_dir, "metadata_barcodes.txt.gz"))
    assignments = meta[meta.columns[-1]].values
    keep = pd.notna(assignments)            # drop cells with no celltype
    assignments = assignments[keep]
    celltypes_sorted = sorted(pd.unique(assignments))
    print(f"{len(celltypes_sorted)} cell types, {keep.sum()} cells", file=sys.stderr)

    # Per-celltype averaged matrix for every method.
    # Optional override: point the sctransform panel at a precomputed residual
    # matrix (e.g. R sctransform v2 float32 .bin) via SCT_OVERRIDE=<path.bin>.
    sct_override = os.environ.get("SCT_OVERRIDE")
    gb = {}
    for label, fname, dense in PANELS:
        if fname == "__CLR__":
            from norm_sparse import norm_clr
            from metrics_matrix import compute_overdispersion
            raw = mmread(os.path.join(dataset_dir, "subset_genes", "raw.mtx.gz")).tocsr()
            alpha = float(compute_overdispersion(raw))
            scale = float(np.asarray(raw.sum(1)).ravel().mean())
            print(f"computing {label} (additive CLR, delta-method K=4*alpha*s, "
                  f"alpha={alpha:.3g})...", file=sys.stderr)
            X = np.asarray(norm_clr(raw, alpha=alpha, scale=scale), dtype=np.float32)[keep]
        else:
            if label == "sctransform" and sct_override:
                path = sct_override
            else:
                path = os.path.join(dataset_dir, "subset_genes", fname)
            print(f"loading {label} from {path}...", file=sys.stderr)
            X = load_matrix(path, dense)[keep]
        gb[label] = celltype_means(X, celltypes_sorted, assignments)
        del X

    # Top-n_top genes per celltype, defined by the reference method (PF).
    ridx = np.argsort(gb[REF], axis=1)[:, -n_top:]

    fig, axs = plt.subplots(figsize=(10, 10), ncols=2, nrows=2)
    im = None
    for (label, _, _), ax in zip(PANELS, axs.ravel()):
        mat = np.take_along_axis(gb[label], ridx, axis=1)
        mat = stats.rankdata(mat, axis=1, method="ordinal")
        im = ax.imshow(mat, aspect="auto", cmap="cividis")
        title = "sctransform v2" if label == "sctransform" else label
        ax.set(ylabel="Celltype", xlabel=f"Top {n_top} genes",
               xticks=[], yticks=[], title=title)

    fig.colorbar(im, ax=axs.ravel().tolist(), shrink=1, label="Gene rank")

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    fig.savefig(f"{out_prefix}.pdf", facecolor="white", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}.png", facecolor="white", dpi=300, bbox_inches="tight")
    print(f"wrote {out_prefix}.pdf + .png")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    main(sys.argv[1], sys.argv[2], n)
