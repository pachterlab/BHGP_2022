#!/usr/bin/env python3
"""Compute shifted CLR (centered log-ratio) normalization on a raw sparse matrix.

Writes a single dense CSV.gz output: <out_prefix>/clr.csv.gz

Mathematically: T(x)_i = log(x_i/s + c) - mean_j(log(x_j/s + c)) with c=1 by default,
where s = sum_j x_j is the per-cell depth. Equivalent to PFlog1pPF where the second
PF step is replaced by additive cell-centering (see docs/current_paper/supplementary-note).

Output is dense because per-cell mean subtraction destroys sparsity.

Usage:
    python norm_clr.py <raw.mtx.gz> <out_prefix_dir> [c]
"""
import os
import sys

import pandas as pd
from scipy.io import mmread

from norm_sparse import norm_clr


def main(in_matrix_fn, out_prefix, c=1.0):
    print(f"loading {in_matrix_fn}", flush=True)
    mtx = mmread(in_matrix_fn).tocsr()
    print(f"  shape: {mtx.shape}, density: {mtx.nnz / (mtx.shape[0]*mtx.shape[1]):.4f}", flush=True)

    print(f"computing CLR (c={c})", flush=True)
    out = norm_clr(mtx, c=c)

    out_fn = os.path.join(out_prefix, "clr.csv.gz")
    print(f"writing {out_fn}", flush=True)
    pd.DataFrame(out).to_csv(out_fn, index=False, header=False, compression="gzip")
    print(f"  done: shape {out.shape}, value range [{out.min():.4f}, {out.max():.4f}]", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    in_fn = args[0]
    out_dir = args[1]
    c = float(args[2]) if len(args) > 2 else 1.0
    main(in_fn, out_dir, c=c)
