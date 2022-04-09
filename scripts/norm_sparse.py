#!/usr/bin/env python3

import sys
import os
from scipy.io import mmread, mmwrite

import numpy as np
import scipy as sp

def do_pf(mtx, sf = None):
    pf = mtx.sum(axis=1).A.ravel()
    if not sf:
        sf = pf.mean()
    pf = sp.sparse.diags(sf/pf) @ mtx
    return pf

def norm_raw(mtx):
    return mtx

def norm_pf(mtx):
    return do_pf(mtx)

def norm_log(mtx):
    return np.log1p(mtx)

def norm_pf_log(mtx):
    pf_log = np.log1p(do_pf(mtx))
    return pf_log

def norm_pf_log_pf(mtx):
    pf_log_pf = do_pf(np.log1p(do_pf(mtx)))
    return pf_log_pf

def norm_cpm_log(mtx):
    cpm_log = np.log1p(do_pf(mtx, sf=1e6))
    return cpm_log

def norm_cp10k_log(mtx):
    cp10k_log =  np.log1p(do_pf(mtx, sf=1e4))
    return cp10k_log

def norm_sqrt(mtx):
    sqrt = np.sqrt(mtx)
    return sqrt

NORM = {
    "raw": norm_raw,
    "pf": norm_pf,
    "log": norm_log,
    "pf_log": norm_pf_log,
    "pf_log_pf": norm_pf_log_pf,
    "cpm_log": norm_cpm_log,
    "cp10k_log": norm_cp10k_log,
    "sqrt": norm_sqrt,
}



def main(in_matrix_fn, out_prefix):
    mtx = mmread(in_matrix_fn)

    titles = ["pf", "log", "pf_log", "pf_log_pf", "cpm_log", "cp10k_log", "sqrt"]
    for title in titles:
        out_fn = os.path.join(out_prefix, f"{title}.mtx")
        mmwrite(out_fn, NORM[title](mtx))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
