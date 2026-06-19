#!/usr/bin/env python3
"""Sample-check whether switching sctransform pysctransform 0.1.1 -> R v2 shifts
the summary-boxplot metrics enough to warrant a full re-run.

For ~N datasets spanning the avg_per_cell range, compute cov_gene / r2_depth /
r_mono for sctransform under BOTH implementations (identical metric code) and
write a paired TSV. R v2 residuals are generated per dataset to a temp float32
.bin (via gen_sctransform_v2.R) and deleted immediately to bound disk use.

Usage:
  python sample_check_sct.py <data_root> <out_tsv> [n_datasets] [r_mono_n_cells]
"""
import glob
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.io import mmread
from scipy import stats
from sklearn.linear_model import LinearRegression

HERE = os.path.dirname(os.path.abspath(__file__))
R_ENV = os.environ.copy()
if os.environ.get("SCTRANSFORM_R_LD_LIBRARY_PATH"):
    R_ENV["R_LD_LIBRARY_PATH"] = os.environ["SCTRANSFORM_R_LD_LIBRARY_PATH"]


def load_bin(path):
    n_cells, n_genes = map(int, open(path + ".shape").read().strip().split(","))
    return np.fromfile(path, dtype=np.float32).reshape(n_cells, n_genes)


def load_csv_gz(path):
    import gzip
    with gzip.open(path, "rt") as f:
        return np.loadtxt(f, delimiter=",", dtype=np.float32)


def cov_gene(X):
    gv = np.var(X, axis=0)
    return float(np.sqrt(np.var(gv)) / np.mean(gv))


def r2_depth(X, raw):
    x = np.asarray(raw.sum(1)).ravel().astype(float)
    y = X.sum(1).ravel().astype(float)
    minx, maxx = x.min(), x.max()
    miny, maxy = y.min(), y.max()
    maxy = maxy - miny
    if maxy == 0:
        return 0.0
    xx = ((x - minx) / maxx).reshape(-1, 1)
    yy = (y - miny) / maxy
    return float(LinearRegression().fit(xx, yy).score(xx, yy))


def r_mono(X, raw, n_sub, seed=0):
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    idx = rng.choice(n, n_sub, replace=False) if n > n_sub else np.arange(n)
    rv = np.empty(len(idx))
    for j, i in enumerate(idx):
        rv[j], _ = stats.spearmanr(X[i], np.asarray(raw[i].todense()).ravel())
    return float(np.nanmean(rv))


def pick_datasets(data_root, n):
    cands = []
    for fn in sorted(glob.glob(os.path.join(data_root, "*", "subset_genes",
                                            "*_subset_genes_metrics.json"))):
        ds = os.path.basename(os.path.dirname(os.path.dirname(fn)))
        d = os.path.dirname(fn)
        if not (os.path.exists(os.path.join(d, "raw.mtx.gz")) and
                os.path.exists(os.path.join(d, "sctransform.csv.gz"))):
            continue
        try:
            apc = json.load(open(fn)).get("avg_per_cell")
        except Exception:
            apc = None
        if apc:
            cands.append((apc, ds, d))
    cands.sort()
    if len(cands) <= n:
        return cands
    step = len(cands) / n
    return [cands[int(i * step)] for i in range(n)]


def main(data_root, out_tsv, n_datasets=30, r_mono_n=1500):
    sel = pick_datasets(data_root, int(n_datasets))
    print(f"selected {len(sel)} datasets", file=sys.stderr)
    rows = ["ds\tavg_per_cell\timpl\tcov_gene\tr2_depth\tr_mono"]
    for apc, ds, d in sel:
        raw = mmread(os.path.join(d, "raw.mtx.gz")).tocsr().astype(np.float32)
        print(f"[{ds}] raw {raw.shape}, apc={apc:.0f}", file=sys.stderr)
        # pysctransform (stored)
        try:
            Xp = load_csv_gz(os.path.join(d, "sctransform.csv.gz"))
            rows.append(f"{ds}\t{apc:.1f}\tpysct\t{cov_gene(Xp):.5f}\t{r2_depth(Xp,raw):.5f}\t{r_mono(Xp,raw,r_mono_n):.5f}")
        except Exception as e:
            print(f"  pysct FAIL: {e}", file=sys.stderr)
        # R v2
        with tempfile.TemporaryDirectory() as td:
            binp = os.path.join(td, f"{ds}.bin")
            p = subprocess.run(["Rscript", os.path.join(HERE, "gen_sctransform_v2.R"),
                                os.path.join(d, "raw.mtx.gz"), binp],
                               env=R_ENV, capture_output=True, text=True)
            if not os.path.exists(binp):
                print(f"  Rv2 FAIL: {p.stderr[-500:]}", file=sys.stderr)
                continue
            Xr = load_bin(binp)
            rows.append(f"{ds}\t{apc:.1f}\trv2\t{cov_gene(Xr):.5f}\t{r2_depth(Xr,raw):.5f}\t{r_mono(Xr,raw,r_mono_n):.5f}")
        open(out_tsv, "w").write("\n".join(rows) + "\n")   # checkpoint each ds
        print(f"  done {ds}", file=sys.stderr)
    print(f"wrote {out_tsv}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
         sys.argv[3] if len(sys.argv) > 3 else 30,
         int(sys.argv[4]) if len(sys.argv) > 4 else 1500)
