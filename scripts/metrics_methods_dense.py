#!/usr/bin/env python3
# compute normalization level metrics

import numpy as np
import pandas as pd
import os
import sys
import json

from scipy.io import mmread
from scipy import stats
from sklearn.linear_model import LinearRegression


def metrics_gcov(mtx):
    gvar = np.var(mtx, axis=0)
    
    gcov = np.sqrt(np.var(gvar))/np.mean(gvar)
    return gcov
    
def metrics_ccov(mtx):
    return None

def metrics_depth(mtx, raw):
    x = raw.sum(1).A.ravel()
    y = mtx.sum(1).ravel() 

    minx, maxx = min(x), max(x)
    miny, maxy = min(y), max(y)
    maxy = maxy - miny

    xx = (x - minx)/maxx
    yy = (y - miny)/maxy

    close = np.all(np.allclose(y, y[0]))
    if close:
        # happens when mtx == PF or PFlog1pPF
        # ie, all cell depths same, no correlation
        r2 = 0
        return r2

    reg = LinearRegression().fit(xx.reshape(-1,1), yy)
    r2 = reg.score(xx.reshape(-1,1), yy)
    return r2
    
def metrics_mono(mtx, raw):
    rv = np.ones(mtx.shape[0])
    for i in range(mtx.shape[0]):
        r, p  = stats.spearmanr(mtx[i], raw[i].A.ravel())
        rv[i] = r
    return rv.mean()
    
    
def compute_method_metrics(mtx, raw):
    entry = {
      'cov_gene': None,
      'cov_cell': None,
      'r2_depth': None,
      'r_mono': None
    }
    try:
        entry["cov_gene"] = metrics_gcov(mtx)
    except:
        pass

    try:
        entry["cov_cell"] = metrics_ccov(mtx)
    except:
        pass

    try:
        entry["r2_depth"] = metrics_depth(mtx, raw)
    except:
        pass

    try:
        entry["r_mono"  ] = metrics_mono(mtx, raw)
    except:
        pass
    return entry

def main(raw_mtx_fn, norm_mtx_fn, metrics_fn):
    txlabel = {
      'raw': 'raw',
      'pf': 'PF',
      "sqrt": "sqrt",
      'log': 'log1p',
      'cp10k_log': 'log1pCP10k',
      'cpm_log': 'log1pCPM',
      'cp10k_log_scale': 'scalelog1pCP10k',
      'sctransform': 'sctransform',
      'pf_log': 'log1pPF',
      'pf_log_pf': 'PFlog1pPF'
    }
    method = txlabel[os.path.basename(norm_mtx_fn).split(".")[0]]
    
    raw = mmread(raw_mtx_fn).tocsr()
    mtx = pd.read_csv(norm_mtx_fn, header=None, compression="gzip").values
    entry = compute_method_metrics(mtx, raw)
    
    # read, update, write
    with open(metrics_fn, "r") as jsonFile:
        metrics = json.load(jsonFile)
    metrics[method] = entry
    with open(metrics_fn, "w") as jsonFile:
        json.dump(metrics, jsonFile)

if __name__ == "__main__":
    # raw_mtx_fn, norm_mtx_fn, metrics_fn
    main(sys.argv[1], sys.argv[2], sys.argv[3])
