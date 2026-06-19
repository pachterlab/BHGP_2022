#!/usr/bin/env python3
"""Recompute the summary metrics for 'PFlog (shift. CLR)' using the ADDITIVE
centered log-ratio (norm_clr), replacing the multiplicative pf_log_pf values in
every dataset's *_subset_genes_metrics.json.

norm_clr is pure Python/scipy, so this needs no R and no persisted clr.csv.gz:
each dataset's CLR is computed on the fly from raw.mtx.gz, metrics taken, matrix
discarded. Other methods' entries are left untouched. Resumable via a manifest.

Usage: python rerun_pflog_clr.py <data_root> [n_workers] [manifest]
"""
import glob, json, os, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.io import mmread
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import norm_sparse as NS          # norm_clr (additive shifted CLR)
import compare_sct_impl as C      # cov_gene / r2_depth / r_mono (same defs as metrics_methods)
from metrics_matrix import compute_overdispersion   # dataset overdispersion alpha

LABEL = "PFlog (shift. CLR)"


def find_json(d):
    j = glob.glob(os.path.join(d, "*_subset_genes_metrics.json"))
    return j[0] if j else None


def process(args):
    ds, d = args
    raw_fn, jpath = os.path.join(d, "raw.mtx.gz"), find_json(d)
    try:
        raw = mmread(raw_fn).tocsr().astype(np.float64)
        # Delta-method first-PF constant K = 4*alpha*s, with the dataset (pooled)
        # overdispersion alpha and scale s = mean total UMI per cell. This is the
        # variance-stabilizing pseudocount y0 = 1/(4*alpha); it lowers cov_gene
        # while leaving r2_depth (~0) and r_mono (=1) unchanged (both K-invariant).
        alpha = float(compute_overdispersion(raw))
        scale = float(np.asarray(raw.sum(1)).ravel().mean())
        X = np.asarray(NS.norm_clr(raw, alpha=alpha, scale=scale))   # additive shifted CLR (dense)
        raw32 = raw.astype(np.float32)
        entry = {"cov_gene": C.cov_gene(X), "cov_cell": None,
                 "r2_depth": C.r2_depth(X, raw32), "r_mono": C.r_mono(X, raw32),
                 "clr_alpha": alpha, "clr_scale": scale, "clr_K": 4.0 * alpha * scale}
        dd = json.load(open(jpath))
        if isinstance(dd.get(LABEL), dict):
            dd[LABEL].update(entry)
        else:
            dd[LABEL] = entry
        dd["pflog_impl"] = "additive CLR (norm_clr), delta-method K=4*alpha*s"
        json.dump(dd, open(jpath, "w"))
        return (ds, "OK", f"cov={entry['cov_gene']:.2f} r2d={entry['r2_depth']:.3f} rmono={entry['r_mono']:.3f} K={entry['clr_K']:.0f}")
    except Exception as e:
        return (ds, "FAIL", repr(e)[:200])


def main(data_root, n_workers=8, manifest="/tmp/pflog_clr_done.txt"):
    n_workers = int(n_workers)
    done = {l.split("\t")[0] for l in open(manifest)} if os.path.exists(manifest) else set()
    tasks = []
    for d in sorted(glob.glob(os.path.join(data_root, "*", "subset_genes"))):
        ds = os.path.basename(os.path.dirname(d))
        if ds in done:
            continue
        if os.path.exists(os.path.join(d, "raw.mtx.gz")) and find_json(d):
            tasks.append((ds, d))
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
