#!/usr/bin/env python3
"""Compute corrected runorm/scclr PFlog normalization on a raw sparse matrix.

Writes a single dense CSV.gz output: <out_prefix>/clr.csv.gz

Mathematically, PFlog is the shifted CLR
center(log1p(4 * alpha * x)), computed by runorm/scclr.  This is not the older
depth-targeted CLR helper center(log1p(K * x / s_i)).

Output is dense because per-cell mean subtraction destroys sparsity.

Usage:
    python norm_clr.py <raw.mtx.gz> <out_prefix_dir>
"""
import os
import sys

import pandas as pd
from scipy.io import mmread

from scclr_pflog import normalize_pflog, to_dense


def main(in_matrix_fn, out_prefix):
    print(f"loading {in_matrix_fn}", flush=True)
    mtx = mmread(in_matrix_fn).tocsr()
    print(f"  shape: {mtx.shape}, density: {mtx.nnz / (mtx.shape[0]*mtx.shape[1]):.4f}", flush=True)

    print("computing PFlog with runorm/scclr target='auto'", flush=True)
    out = to_dense(normalize_pflog(mtx))

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
    if len(args) > 2:
        print("WARNING: the legacy c argument is ignored; PFlog uses scclr target='auto'", file=sys.stderr)
    main(in_fn, out_dir)
