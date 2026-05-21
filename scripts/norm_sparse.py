#!/usr/bin/env python3

import sys
import os
from scipy.io import mmread, mmwrite

import numpy as np
import scipy as sp

def do_pf(mtx, sf = None):
    pf = np.asarray(mtx.sum(axis=1)).ravel()
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

def norm_clr(mtx, c=1.0):
    """Shifted CLR: T(x)_i = log(x_i/s + c) - mean_j(log(x_j/s + c)).

    Output is dense — per-cell mean centering destroys sparsity. With c=1, this
    is mathematically equivalent to PFlog1pPF where the second PF is additive
    cell-centering (see docs/current_paper/supplementary-note/supplementary-note.tex
    for the proof). Inherits Aitchison's scale invariance for free; per-cell sums
    of the output are exactly 0 for any c > 0.
    """
    # u = x / s (each row sums to 1.0). Use sf=1.0 explicitly; do_pf's truthy
    # check (`if not sf`) treats 1.0 as truthy so this branch is preserved.
    u = do_pf(mtx, sf=1.0)
    D = mtx.shape[1]
    log_c = np.log(c)

    # Sparse "offset trick": store log(u_nz + c) - log(c) at nonzero positions.
    # That way zero entries (where u=0) hold the natural sparse value 0, and the
    # true log(u+c) at any position is (stored_value + log(c)).
    log_term = u.copy()
    if c == 1.0:
        # log(u + 1) - log(1) = log1p(u); use log1p for numerical stability.
        log_term.data = np.log1p(log_term.data)
    else:
        log_term.data = np.log(log_term.data + c) - log_c

    # Per-cell mean of TRUE values (dense, includes zero positions):
    #   sum_true = sum_stored + nnz*log(c) (true at nonzero)
    #            + (D - nnz)*log(c) (true at zero)
    #            = sum_stored + D*log(c)
    #   mean_true = sum_stored/D + log(c)
    row_sum = np.asarray(log_term.sum(axis=1)).ravel()
    per_cell_mean = row_sum / D + log_c

    # Materialize dense; nonzero stored values need + log(c), zero positions
    # need log(c) (since their true value is log(0+c) = log(c)).
    dense_log = log_term.toarray() + log_c
    return dense_log - per_cell_mean[:, None]

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
