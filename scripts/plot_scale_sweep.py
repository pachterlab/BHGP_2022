#!/usr/bin/env python3
"""Scale-invariance figure: relative Frobenius distance vs total-count scaling.

For total-count scale factors s (default 0.1x .. 10x), scale the raw counts by s,
apply a normalization f, and plot the RELATIVE Frobenius distance to the unscaled
(s=1) baseline:  ||f(sM) - f(M)||_F / ||f(M)||_F.

A normalization that is invariant to sequencing depth / total-count scaling stays
flat near zero; one that is not grows as s moves away from 1.

raw/log1pPF/PFlog use continuous s*M (defined for real input, isolating pure
scaling). sctransform requires integer counts, so its input is round(sM); at
s<1 that rounding also downsamples, so its low-s points reflect information loss
on top of (non-)invariance -- worth a caption note.

Methods:
  - raw                  : no normalization (identity)      -> grows ~ |s-1|
  - log1pPF              : log1p(PF)                         -> grows
  - PFlog (shift. CLR)   : scclr center(log1p(4*alpha*x))
  - sctransform          : Satija v2 Pearson residuals (computed in R via glmGamPoi)

All methods run on the SAME subsampled, sanitized matrix. sctransform is delegated
to scale_sweep_sctransform.R; set SCALE_SWEEP_R_LD_LIBRARY_PATH if the R runtime
needs a custom library path.

Usage:
    python plot_scale_sweep.py <raw.mtx.gz> <out_prefix> [scales_csv] [n_cells] [seed]
"""
import os
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 15
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy import sparse

from scclr_pflog import normalize_pflog, to_dense

R_ENV = os.environ.copy()
if os.environ.get("SCALE_SWEEP_R_LD_LIBRARY_PATH"):
    R_ENV["R_LD_LIBRARY_PATH"] = os.environ["SCALE_SWEEP_R_LD_LIBRARY_PATH"]

cividis = matplotlib.colormaps["cividis"]
CLR_LABEL = "PFlog (shift. CLR)"
STYLE = {
    "raw":         dict(color="k",          marker="o"),
    "log1pPF":     dict(color=cividis(0.5), marker="s"),
    CLR_LABEL:     dict(color="#E41A1C",    marker="D"),
    "sctransform": dict(color=cividis(0.0), marker="^"),
}
ORDER = ["raw", "log1pPF", CLR_LABEL, "sctransform"]


def do_pf(X, sf=None):
    """Proportional fitting to the matrix mean depth (matches norm_sparse.do_pf)."""
    pf = X.sum(axis=1)
    if sf is None:
        sf = pf.mean()
    return X * (sf / pf)[:, None]


def f_raw(X):
    return X


def f_log1pPF(X):
    return np.log1p(do_pf(X))


def f_clr(X):
    """Current scclr PFlog: center(log1p(4*alpha*x))."""
    return to_dense(normalize_pflog(sparse.csr_matrix(X)), dtype=np.float64)


PY_METHODS = {"raw": f_raw, "log1pPF": f_log1pPF, CLR_LABEL: f_clr}


def apply_full(f, scaled):
    """Apply f to the nonzero-sum cells only, returning a full (cells x genes)
    matrix with all-zero rows for cells that scaling+rounding zeroed out. Keeps
    every scale the same shape (matching the R sctransform convention) and avoids
    divide-by-zero in PF on empty cells."""
    out = np.zeros_like(scaled)
    nz = scaled.sum(1) > 0
    if nz.any():
        out[nz] = f(scaled[nz])
    return out


def rel_frob(f1, fs):
    base = np.sqrt((f1 ** 2).sum())
    return float(np.sqrt(((fs - f1) ** 2).sum()) / base) if base > 0 else np.nan


def main(raw_fn, out_prefix, scales_csv="0.1,0.2,0.5,1,2,5,10", n_cells=3000, seed=0):
    scales = [float(x) for x in scales_csv.split(",")]
    rng = np.random.default_rng(seed)

    m = mmread(raw_fn)
    M = np.asarray(m.todense(), dtype=np.float64) if sparse.issparse(m) else np.asarray(m, np.float64)
    # cells x genes. Subsample cells, then sanitize (drop all-zero genes/cells).
    if M.shape[0] > n_cells:
        M = M[np.sort(rng.choice(M.shape[0], n_cells, replace=False))]
    M = M[:, M.sum(0) > 0]
    M = M[M.sum(1) > 0]
    print(f"matrix: {M.shape[0]} cells x {M.shape[1]} genes", file=sys.stderr)

    # --- Python methods (relative Frobenius vs s=1) ---
    # Continuous scaling (no rounding): these transforms are defined for real-
    # valued input, so s*M isolates pure total-count scaling. Only sctransform
    # (below) rounds, since it requires integer counts.
    curves = {k: [] for k in PY_METHODS}
    base = {k: apply_full(f, M) for k, f in PY_METHODS.items()}
    for s in scales:
        scaled = M * s
        for k, f in PY_METHODS.items():
            curves[k].append(rel_frob(base[k], apply_full(f, scaled)))
        print(f"  s={s:g} " + " ".join(f"{k}={curves[k][-1]:.4f}" for k in PY_METHODS),
              file=sys.stderr)

    # --- sctransform v2 via R, on the identical sanitized matrix ---
    sct = None
    with tempfile.TemporaryDirectory() as td:
        raw_tmp = os.path.join(td, "M0.mtx")
        mmwrite(raw_tmp, sparse.csr_matrix(M))   # cells x genes; R transposes
        sct_tsv = os.path.join(td, "sct.tsv")
        r_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "scale_sweep_sctransform.R")
        cmd = ["Rscript", r_script, raw_tmp, scales_csv, sct_tsv,
               str(M.shape[0] + 1), str(seed)]   # n_cells > rows -> no resubsample
        print("running sctransform sweep in R ...", file=sys.stderr)
        p = subprocess.run(cmd, env=R_ENV, capture_output=True, text=True)
        if p.returncode != 0 or not os.path.exists(sct_tsv):
            print("WARNING: sctransform R sweep failed; plotting without it.", file=sys.stderr)
            print(p.stderr[-2000:], file=sys.stderr)
        else:
            sct = pd.read_csv(sct_tsv, sep="\t")

    # --- plot ---
    fig, ax = plt.subplots(figsize=(6, 5))
    for k in ORDER:
        if k in curves:
            ax.plot(scales, curves[k], label=k, markersize=6, **STYLE[k])
        elif k == "sctransform" and sct is not None:
            ax.plot(sct["scale"], sct["rel_frob"], label="sctransform (v2)",
                    markersize=6, **STYLE[k])
    ax.axvline(1.0, color="grey", lw=0.8, ls=":")
    ax.set(xscale="log", xlabel="Total-count scale factor",
           ylabel=r"Relative Frobenius distance")
    ax.legend(frameon=False, fontsize=12)

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    fig.savefig(f"{out_prefix}.pdf", facecolor="white", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}.png", facecolor="white", dpi=300, bbox_inches="tight")
    print(f"wrote {out_prefix}.pdf + .png")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    a = sys.argv
    main(a[1], a[2],
         a[3] if len(a) > 3 else "0.1,0.2,0.5,1,2,5,10",
         int(a[4]) if len(a) > 4 else 3000,
         int(a[5]) if len(a) > 5 else 0)
