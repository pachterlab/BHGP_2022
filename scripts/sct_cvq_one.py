#!/usr/bin/env python3
"""Per-dataset: empirical CV(q) for sctransform from the FULL-vst stored residuals
(sctransform.csv.gz = full-fit Pearson residuals), computed on the filtered genes
(mu>0.05 & detected>=1%) to match the functional CV(q) panel.

  cv_all  = CV across all genes of per-gene residual variance  (== stored cov_gene; sanity)
  cv_filt = CV across filtered genes of per-gene residual variance  (panel-comparable)

Usage: sct_cvq_one.py <dataset> <outdir>
Writes one CSV line to <outdir>/<dataset>.csv (skips if it already exists)."""
import sys, os, json, numpy as np, pandas as pd, scipy.io as sio, warnings
warnings.filterwarnings("ignore")
D = os.environ.get("BHGP_DATA_ROOT", "data")
ds, outdir = sys.argv[1], sys.argv[2]
out = os.path.join(outdir, ds + ".csv")
if os.path.exists(out):
    sys.exit(0)
try:
    m = json.load(open(f"{D}/{ds}/subset_genes/{ds}_subset_genes_metrics.json"))
    nc, ng, a = m["ncells"], m["ngenes"], float(m["overdispersion"])
    stored = float(m["sctransform"]["cov_gene"])
    R = sio.mmread(f"{D}/{ds}/subset_genes/raw.mtx.gz").tocsr().astype(float)
    if R.shape[0] == ng and R.shape[1] == nc:
        R = R.T.tocsr()                                   # cells x genes
    mu = np.asarray(R.mean(0)).ravel(); det = np.asarray((R > 0).sum(0)).ravel() / nc
    mask = (mu > 0.05) & (det >= 0.01)
    F = f"{D}/{ds}/subset_genes/sctransform.csv.gz"
    s = s2 = None; n = 0
    for ck in pd.read_csv(F, header=None, dtype=np.float32, chunksize=2000):
        X = ck.values.astype(np.float64)                  # rows=cells, cols=genes
        if s is None: s = np.zeros(X.shape[1]); s2 = np.zeros(X.shape[1])
        s += X.sum(0); s2 += (X * X).sum(0); n += X.shape[0]
    if s.shape[0] != ng:                                  # orientation/gene-count guard
        raise ValueError(f"residual cols {s.shape[0]} != ngenes {ng}")
    var = s2 / n - (s / n) ** 2
    def CV(v): v = v[np.isfinite(v) & (v > 0)]; return float(np.sqrt(np.var(v)) / np.mean(v))
    cv_all, cv_filt = CV(var), CV(var[mask])
    with open(out, "w") as g:
        g.write(f"{ds},{nc},{ng},{int(mask.sum())},{a:.5f},{cv_all:.5f},{cv_filt:.5f},{stored:.5f}\n")
except Exception as e:
    with open(out + ".err", "w") as g:
        g.write(f"{ds},{type(e).__name__}: {e}\n")
