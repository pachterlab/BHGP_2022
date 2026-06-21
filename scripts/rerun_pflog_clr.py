#!/usr/bin/env python3
"""Recompute only the summary metrics for 'PFlog (shift. CLR)' with scclr.

PFlog is the shifted CLR ``center(log1p(4 * alpha * x))``.  Each
dataset's PFlog object is computed on the fly from raw.mtx.gz, summary metrics
are updated in-place, and other methods' entries are left untouched.  Resumable
via a manifest.

Usage: python rerun_pflog_clr.py <data_root> [n_workers] [manifest] [scclr_path]
"""
import glob, json, os, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.io import mmread
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from metrics_matrix import compute_overdispersion   # dataset overdispersion alpha
from scclr_pflog import cov_gene, normalize_pflog

LABEL = "PFlog (shift. CLR)"
OLD_LABEL = "PFlogPF (shift. CLR)"


def find_json(d):
    j = glob.glob(os.path.join(d, "*_subset_genes_metrics.json"))
    return j[0] if j else None


def process(args):
    ds, d, scclr_path = args
    raw_fn, jpath = os.path.join(d, "raw.mtx.gz"), find_json(d)
    try:
        raw = mmread(raw_fn).tocsr().astype(np.float64)
        alpha = float(compute_overdispersion(raw))
        sclr = normalize_pflog(raw, alpha=alpha, scclr_path=scclr_path)
        entry = {
            "cov_gene": cov_gene(sclr),
            "cov_cell": None,
            "r2_depth": 0.0,
            "r_mono": 1.0,
            "scclr_alpha": float(sclr.alpha),
            "scclr_k": float(sclr.k),
            "pflog_formula": "center(log1p(4*alpha*x))",
        }
        dd = json.load(open(jpath))
        dd.pop(OLD_LABEL, None)
        dd.pop("pflogpf_impl", None)
        dd[LABEL] = entry
        dd["pflog_impl"] = "scclr PFlog: center(log1p(4*alpha*x))"
        json.dump(dd, open(jpath, "w"))
        return (ds, "OK", f"cov={entry['cov_gene']:.2f} r2d=0.000 rmono=1.000 alpha={entry['scclr_alpha']:.3g}")
    except Exception as e:
        return (ds, "FAIL", repr(e)[:200])


def main(data_root, n_workers=8, manifest="/tmp/pflog_scclr_done.txt", scclr_path=None):
    n_workers = int(n_workers)
    done = {l.split("\t")[0] for l in open(manifest)} if os.path.exists(manifest) else set()
    tasks = []
    for d in sorted(glob.glob(os.path.join(data_root, "*", "subset_genes"))):
        ds = os.path.basename(os.path.dirname(d))
        if ds in done:
            continue
        if os.path.exists(os.path.join(d, "raw.mtx.gz")) and find_json(d):
            tasks.append((ds, d, scclr_path))
    print(f"{len(tasks)} datasets ({len(done)} done), {n_workers} workers", flush=True)
    ok = fail = 0
    with ProcessPoolExecutor(max_workers=n_workers) as ex, open(manifest, "a") as mf:
        futs = {ex.submit(process, t): t[0] for t in tasks}
        for i, fut in enumerate(as_completed(futs), 1):
            ds, st, info = fut.result()
            ok += st == "OK"; fail += st != "OK"
            mf.write(f"{ds}\t{st}\t{info}\n"); mf.flush()
            print(f"[{i}/{len(tasks)}] {ds} {st} {info}", flush=True)
    print(f"DONE: {ok} ok, {fail} failed", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
