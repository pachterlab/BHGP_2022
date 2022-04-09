#!/usr/bin/env python3
# compute normalization level metrics

import numpy as np
import os
import sys
import json

from scipy.io import mmread
from sklearn.linear_model import LinearRegression


def metrics_gcov(mtx):
    mtx_squared = mtx.copy()
    mtx_squared.data **= 2
    gvar =  (mtx_squared.mean(0) - np.square(mtx.mean(0))).A.ravel()
    
    gcov = np.sqrt(np.var(gvar))/np.mean(gvar)
    return gcov
    
def metrics_ccov(mtx):
    csum = mtx.sum(1).A.ravel()
    
    ccov = np.sqrt(np.var(csum))/np.mean(csum)
    return ccov

def metrics_depth(mtx, raw, same=False):
    x = raw.sum(1).A.ravel()
    y = mtx.sum(1).A.ravel() 

    minx, maxx = min(x), max(x)
    miny, maxy = min(y), max(y)
    maxy = maxy - miny

    xx = (x - minx)/maxx
    yy = (y - miny)/maxy
    if same:
        r2 = 1
        # happens when mtx == raw
        return r2

    close = np.all(np.allclose(y, y[0]))
    if close:
        # happens when mtx == PF or PFlog1pPF
        # ie, all cell depths same, no correlation
        r2 = 0
        return r2

    reg = LinearRegression().fit(xx.reshape(-1,1), yy)
    r2 = reg.score(xx.reshape(-1,1), yy)
    return r2
    
def metrics_mono(mtx):
    return 1
    
def compute_method_metrics(mtx, raw, method):
    same = False
    if method == "raw":
        same = True
    entry = {
      'cov_gene': None,
      'cov_cell': None,
      'r2_depth': None,
      'r_mono': None
    }
    try:
        entry["cov_gene"] = metrics_gcov(mtx)
        entry["cov_cell"] = metrics_ccov(mtx)
        entry["r2_depth"] = metrics_depth(mtx, raw, same)
        entry["r_mono"  ] = metrics_mono(mtx)
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
    
    raw = mmread(raw_mtx_fn)
    mtx = mmread(norm_mtx_fn)
    entry = compute_method_metrics(mtx, raw, method)
    
    # read, update, write
    with open(metrics_fn, "r") as jsonFile:
        metrics = json.load(jsonFile)
    metrics[method] = entry
    with open(metrics_fn, "w") as jsonFile:
        json.dump(metrics, jsonFile)

if __name__ == "__main__":
    # raw_mtx_fn, norm_mtx_fn, metrics_fn
    main(sys.argv[1], sys.argv[2], sys.argv[3])
