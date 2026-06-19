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

def norm_clr(mtx, c=1.0, sf=None, alpha=None, scale=None):
    """PFlog (shifted CLR): first PF to depth `sf`, then log(.+c), then additive
    per-cell centering (the CLR step).

    Default sf=None -> mean cell depth, i.e. PF -> log1p -> center. This keeps the
    log compression (variance stabilization) AND the centering (depth removal +
    within-cell rank preservation); output is dense (centering fills zeros) with
    per-cell sums exactly 0. It is a shifted CLR on the composition with effective
    shift c/mean_depth -- a small, variance-stabilization-motivated shift (cf. the
    Supplementary Note's shift-selection criterion).

    If `alpha` (the dataset overdispersion) is given, the first-PF target K is set
    by the delta-method rule K = 4*alpha*scale, where `scale` is the dataset scale
    factor (mean total UMI per cell; defaults to the mean cell depth of `mtx`).
    This calibrates the count-scale pseudocount to y0 = 1/(4*alpha) -- the
    variance-stabilizing pseudocount of the Supplementary Note -- and OVERRIDES sf.
    Because K only rescales the argument of the log (a per-cell monotone map) and
    the additive centering is unchanged, this moves only the variance-stabilization
    behavior; depth removal and within-cell ranks are unaffected.

    sf=1.0 (each cell sums to 1) gives an *exactly* scale-invariant CLR, but with
    c=1 the pseudocount dwarfs the composition (u_i ~ 1e-3 so log(u+1)~u), giving
    no variance stabilization; sf=mean restores it at the cost of only approximate
    global-scale invariance (depth removal via centering is still exact).
    """
    if alpha is not None:
        if scale is None:
            scale = np.asarray(mtx.sum(axis=1)).ravel().mean()
        sf = 4.0 * alpha * scale       # K = 4*alpha*s  (delta-method PF constant)
    # sf=None -> do_pf uses the mean cell depth.
    u = do_pf(mtx, sf=sf)
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
