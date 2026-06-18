#!/usr/bin/env python3
"""Compute CV(q) per functional-form transform for every dataset, on the subset_genes
gene set with the low-info filter (mu>0.05 AND detected in >=1% of cells).

q_g = h'(mu_g)^2 * (mu_g + alpha*mu_g^2)   (model-based, biology-free delta-method
technical variance), CV(q) = sd_g(q)/mean_g(q) across the filtered genes.

scalelog1pCP10k is not a pure abundance map; its q uses the delta-method formula
  q_scale,g = h'_CP10k(mu_g)^2 (mu_g+alpha mu_g^2) / Var_emp(log1pCP10k gene g).

sctransform is omitted: by construction it standardizes the technical variance, so
its CV(q) is 0 (plotted as a reference line, not a functional transform).

Output: a tidy wide CSV with one row per dataset, one CV(q) column per method.
"""
import json, os, glob, gc, sys, numpy as np, scipy.io as sio, scipy.sparse as sp, warnings
warnings.filterwarnings("ignore")

DATA = "/home/sina/projects/synchromesh/data"
OUT  = "/home/sina/projects/BHGP_2022/analysis/figures/cvq_per_method_subset_filtered.csv"
METHODS = ["raw","PF","sqrt","log1p","log1pCP10k","log1pCPM","scalelog1pCP10k","log1pPF","PFlogPF"]

def CVq(a):
    a = a[np.isfinite(a) & (a > 0)]
    return float(np.sqrt(np.var(a)) / np.mean(a)) if a.size else float("nan")

def load(ds):
    mj = f"{DATA}/{ds}/subset_genes/{ds}_subset_genes_metrics.json"
    m = json.load(open(mj)); nc, ng = m["ncells"], m["ngenes"]; a = float(m["overdispersion"])
    X = sio.mmread(f"{DATA}/{ds}/subset_genes/raw.mtx.gz").tocsr().astype(float)
    if X.shape[0] == ng and X.shape[1] == nc:
        X = X.T.tocsr()                                   # -> cells x genes
    return X, a

def cvq_row(X, a):
    nc = X.shape[0]
    s = np.asarray(X.sum(1)).ravel(); sbar = float(s.mean())
    mu = np.asarray(X.mean(0)).ravel()
    det = np.asarray((X > 0).sum(0)).ravel() / nc
    mask = (mu > 0.05) & (det >= 0.01)
    muf = mu[mask]; V = muf + a * muf**2
    g10, gM = 1e4 / sbar, 1e6 / sbar
    # empirical per-gene variance of log1p(CP10k) for scale (sparse: zeros stay zero)
    L = sp.diags(1e4 / s).dot(X).tocsr(); L.data = np.log1p(L.data)
    EL = np.asarray(L.sum(0)).ravel() / nc
    EL2 = np.asarray(L.multiply(L).sum(0)).ravel() / nc
    var_emp = (EL2 - EL**2)[mask]
    hc = g10 / (1 + g10 * muf)
    H = {"raw": np.ones_like(muf), "PF": np.ones_like(muf), "sqrt": 0.5/np.sqrt(muf),
         "log1p": 1/(1+muf), "log1pCP10k": hc, "log1pCPM": gM/(1+gM*muf),
         "log1pPF": 1/(1+muf), "PFlogPF": (4*a)/(1+4*a*muf)}
    row = {k: CVq(H[k]**2 * V) for k in H}
    row["scalelog1pCP10k"] = CVq((hc**2 * V) / var_emp)
    return row, int(mask.sum()), sbar

datasets = sorted(os.path.basename(p).replace("_subset_genes_metrics.json","")
                  for p in glob.glob(f"{DATA}/*/subset_genes/*_subset_genes_metrics.json"))
datasets = [d for d in datasets if os.path.exists(f"{DATA}/{d}/subset_genes/raw.mtx.gz")]
print(f"{len(datasets)} datasets with subset_genes/raw.mtx.gz", flush=True)

out = open(OUT, "w")
out.write("ds,alpha,ncells,ngenes,ngenes_filt,sbar," + ",".join(METHODS) + "\n")
ok = fail = 0
for i, ds in enumerate(datasets):
    try:
        X, a = load(ds)
        if not (np.isfinite(a) and a > 0):
            print(f"  skip {ds}: bad alpha {a}", flush=True); fail += 1; continue
        row, nf, sbar = cvq_row(X, a)
        out.write(f"{ds},{a:.5f},{X.shape[0]},{X.shape[1]},{nf},{sbar:.3f}," +
                  ",".join(f"{row[m]:.5f}" for m in METHODS) + "\n"); out.flush()
        ok += 1
        del X; gc.collect()
    except Exception as e:
        print(f"  FAIL {ds}: {type(e).__name__}: {e}", flush=True); fail += 1
    if (i+1) % 25 == 0:
        print(f"  {i+1}/{len(datasets)}  (ok={ok} fail={fail})", flush=True)
out.close()
print(f"done -> {OUT}  (ok={ok} fail={fail})", flush=True)
