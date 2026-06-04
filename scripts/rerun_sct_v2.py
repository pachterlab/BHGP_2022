#!/usr/bin/env python3
"""Migrate all datasets' sctransform to the R reference v2 (glmGamPoi) and
refresh their summary metrics.

For every dataset under <data_root>/*/subset_genes with raw.mtx.gz + a metrics
json, this:
  1. runs gen_sctransform_v2.R to (over)write subset_genes/sctransform.csv.gz
     with R v2 residuals (the old pysctransform file is moved to
     sctransform.pysct.csv.gz.bak the first time, for reversibility),
  2. recomputes cov_gene / r2_depth / r_mono on the new residuals (same code as
     metrics_methods_dense.py) and updates the 'sctransform' entry in the
     dataset's *_subset_genes_metrics.json (tagging sctransform_impl),
  3. records the dataset in a done-manifest so the run is resumable.

Parallel across datasets; each R process is capped to a few threads.

Usage:
  python rerun_sct_v2.py <data_root> [n_workers] [r_threads] [manifest]
"""
import glob
import gzip
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import stats
from sklearn.linear_model import LinearRegression

HERE = os.path.dirname(os.path.abspath(__file__))
BUILDTOOLS = "/home/sina/bin/miniconda3/envs/buildtools"
GEN_R = os.path.join(HERE, "gen_sctransform_v2.R")
IMPL_TAG = "R sctransform v2 (glmGamPoi)"


def _env(r_threads):
    e = dict(os.environ, R_LD_LIBRARY_PATH=f"/usr/lib/R/lib:{BUILDTOOLS}/lib")
    e["OMP_NUM_THREADS"] = str(r_threads)
    return e


def load_csv_gz(path):
    # pandas C parser is ~10x faster than np.loadtxt for these large files.
    return pd.read_csv(path, header=None, dtype=np.float32).values


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
    return float(LinearRegression().fit(xx, (y - miny) / maxy).score(xx, (y - miny) / maxy))


def r_mono(X, raw):
    rv = np.empty(X.shape[0])
    for i in range(X.shape[0]):
        rv[i], _ = stats.spearmanr(X[i], np.asarray(raw[i].todense()).ravel())
    return float(np.nanmean(rv))


def find_json(d):
    j = glob.glob(os.path.join(d, "*_subset_genes_metrics.json"))
    return j[0] if j else None


def process(args):
    ds, d, r_threads = args
    raw_fn = os.path.join(d, "raw.mtx.gz")
    sct_fn = os.path.join(d, "sctransform.csv.gz")
    jpath = find_json(d)
    try:
        # 1. back up pysct once, then regenerate sctransform.csv.gz with R v2
        bak = os.path.join(d, "sctransform.pysct.csv.gz.bak")
        if os.path.exists(sct_fn) and not os.path.exists(bak):
            shutil.move(sct_fn, bak)
        p = subprocess.run(["Rscript", GEN_R, raw_fn, sct_fn, str(r_threads)],
                           env=_env(r_threads), capture_output=True, text=True)
        if not os.path.exists(sct_fn):
            return (ds, "FAIL_R", p.stderr[-400:])
        # 2. metrics on the new residuals
        X = load_csv_gz(sct_fn)
        raw = mmread(raw_fn).tocsr().astype(np.float32)
        if not np.isfinite(X).all():
            # Degenerate sctransform fit (non-finite residuals) -- same family
            # of pathologically-sparse datasets pysctransform produced all-null
            # for. Record nulls so the summary's NaN-drop excludes them, as before.
            entry = {"cov_gene": None, "cov_cell": None, "r2_depth": None, "r_mono": None}
            status_info = "non-finite residuals -> null metrics"
        else:
            entry = {"cov_gene": cov_gene(X), "cov_cell": None,
                     "r2_depth": r2_depth(X, raw), "r_mono": r_mono(X, raw)}
            status_info = None
        if jpath:
            dd = json.load(open(jpath))
            if isinstance(dd.get("sctransform"), dict):
                dd["sctransform"].update(entry)
            else:
                dd["sctransform"] = entry
            dd["sctransform_impl"] = IMPL_TAG
            json.dump(dd, open(jpath, "w"))
        if status_info:
            return (ds, "OK_NULL", status_info)
        return (ds, "OK", f"cov={entry['cov_gene']:.3f} r2d={entry['r2_depth']:.3f} rmono={entry['r_mono']:.3f}")
    except Exception as e:
        return (ds, "FAIL_PY", repr(e)[:300])


def main(data_root, n_workers=12, r_threads=2, manifest="/tmp/sct_rerun_done.txt"):
    n_workers, r_threads = int(n_workers), int(r_threads)
    done = set()
    if os.path.exists(manifest):
        done = {l.split("\t")[0] for l in open(manifest) if l.strip()}
    tasks = []
    for d in sorted(glob.glob(os.path.join(data_root, "*", "subset_genes"))):
        ds = os.path.basename(os.path.dirname(d))
        if ds in done:
            continue
        if os.path.exists(os.path.join(d, "raw.mtx.gz")) and find_json(d):
            tasks.append((ds, d, r_threads))
    print(f"{len(tasks)} datasets to process ({len(done)} already done), "
          f"{n_workers} workers x {r_threads} threads", flush=True)

    n_ok = n_fail = 0
    with ProcessPoolExecutor(max_workers=n_workers) as ex, open(manifest, "a") as mf:
        futs = {ex.submit(process, t): t[0] for t in tasks}
        for i, fut in enumerate(as_completed(futs), 1):
            ds, status, info = fut.result()
            if status == "OK":
                n_ok += 1
            else:
                n_fail += 1
            mf.write(f"{ds}\t{status}\t{info}\n"); mf.flush()
            print(f"[{i}/{len(tasks)}] {ds} {status} {info}", flush=True)
    print(f"DONE: {n_ok} ok, {n_fail} failed", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
