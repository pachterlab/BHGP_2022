#!/usr/bin/env python3

import os
import sys
import json

from collections import defaultdict

import pandas as pd
import numpy as np

from sklearn.preprocessing import normalize, scale
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from collections import OrderedDict

# from synchromesh.scripts.utils import read_str_list, sanitize_mtx, norm, do_pf, do_log_pf
# from synchromesh.scripts.plot import  plot_depth_norm, plot_depth_dist, plot_knee, plot_pc_depth, plot_mean_var, plot_monotone, plot_example_gene

from scipy.sparse import csr_matrix
from scipy.io import mmread, mmwrite
from scipy import stats
from scipy.sparse import issparse

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import matplotlib.gridspec as gridspec

def nd(arr):
    return np.asarray(arr).reshape(-1)

def yex(ax):
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    
    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    return ax

fsize=15


alpha = 0.33

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams["font.size"] = fsize
#%config InlineBackend.figure_format = 'retina'

cividis = matplotlib.colormaps["cividis"]
colors = {
    "cell": cividis(0.01),
    "gene": cividis(0.5),
    "mono": cividis(0.99)
}

def myvar(a, sparse=True, axis=None):
    """ Variance of sparse matrix a
    var = mean(a**2) - mean(a)**2
    """
    if sparse:
        a_squared = a.copy()
        a_squared.data **= 2
        return np.asarray(a_squared.mean(axis) - np.square(a.mean(axis))).ravel()
    else:
        return np.var(a, axis=axis)
    
def mymean(a, sparse=True, axis=None):
    if sparse:
        return np.asarray(a.mean(axis=axis)).ravel()
    else:
        return a.mean(axis=axis)
    
def mysum(a, sparse=True, axis=None):
    if sparse:
        return np.asarray(a.sum(axis=axis)).ravel()
    else:
        return a.sum(axis=axis)
    
def mygetidx(a, idx, sparse=True, axis=0):
    if sparse:
        if axis == 0:
            return a.getrow(idx).toarray().ravel()
        elif axis == 1:
            return a.getcol(idx).toarray().ravel()
    else:
        if axis == 0:
            return a[idx]
        elif axis == 1:
            return a[idx].ravel()

def get_min_diff(matrix):
    x = matrix.flatten()
    xs = np.sort(x)
    b = xs[0]
    mn = 1100.
    for i in range(1, xs.shape[0]):
        a = xs[i]
        if a == b:
            continue

        diff = abs(a-b)
        if diff < mn and diff > 0:
            mn = diff
        b = a
    return mn

def plot_meanvar(mtx, raw_mean, minlim = 1e-4, maxlim = 1e5, ax=None, gvar=None, cq=None):
    p = {
        "xlabel": "Gene mean",
        "ylabel": "Gene variance",
        "xscale": "log",
        "yscale": "log",
        "xlim": (minlim, maxlim),
    }

    # gvar may be precomputed (e.g. PFlog from the sparse-plus-rank-one form) to
    # avoid densifying the transformed matrix.
    if gvar is None:
        gvar = myvar(mtx, sparse=issparse(mtx), axis=0)
    gcov = np.sqrt(np.var(gvar))/np.mean(gvar)
    # The mean-variance scatter still shows the empirical gene variances, but the
    # reported number is CV(q) (the delta-method metric, when available); the label
    # prefix stays "CoV:" for continuity with the rest of the supplement.
    val = gcov if cq is None else cq

    y = gvar
    yy = (y-y.mean())/np.sqrt(np.var(y))

    ax.scatter(raw_mean, y, facecolor=colors["gene"], alpha=alpha, edgecolor="k", label=f"CoV: {val:,.2f}")
    ax.legend(prop={"size": 12})
    ax.set(**p)
    yex(ax)
    return (ax, val)

def plot_depth(mtx, raw_cell_counts, ax, cell_sums=None):
    x = raw_cell_counts
    # cell_sums may be precomputed (PFlog is centered, so per-cell sums are 0).
    y = mysum(mtx, sparse=issparse(mtx), axis=1) if cell_sums is None else cell_sums

    minx, maxx = min(x), max(x)
    miny, maxy = min(y), max(y)
    maxy = maxy - miny

    xx = (x - minx)/maxx

    # Degenerate case (PF / PFlog): all transformed cell sums equal -> no depth.
    # Handle before the y-rescale to avoid a 0/0 divide.
    close = np.all(np.allclose(y, y[0]))
    if close:
        yy = np.ones(len(y))
    else:
        yy = (y - miny)/maxy
    ax.scatter(xx,yy, edgecolor="k", facecolor=colors["cell"], alpha=alpha)
        
    reg = LinearRegression().fit(xx.reshape(-1,1), yy)
    r2 = reg.score(xx.reshape(-1,1), yy)

    if close:
        # handle the degenerate case where the slope is 0 since all values y are same
        r2 = 0
    
    xxx = np.array([min(xx), max(xx)])

    ax.plot(xxx, reg.coef_*xxx+ reg.intercept_, color="darkgray", linestyle="--", label=f"r$^2$: {r2:,.2f}", linewidth=3)
    

    p = {
      "xlabel": "Raw cell count",
      "ylabel": "Transform cell count",
      "xlim": (-0.1, 1.1),
      "ylim": (-0.1, 1.1),
    }
    ax.set(**p)
    ax.legend(prop={"size": 12})
    return (ax, r2)

def mono(matrix, raw, sample=5000, seed=0):
    sparse = issparse(matrix)
    n = matrix.shape[0]
    if sparse:
        # sparse-preserving transforms are monotonic by construction; per-cell
        # Spearman vs raw is always 1.0 in that case.
        return np.ones(n)
    # Dense transforms need the per-cell Spearman loop, which is the figure's
    # bottleneck on large datasets. A representative cell subsample gives the same
    # monotonicity histogram/mean.
    cells = (np.sort(np.random.default_rng(seed).choice(n, sample, replace=False))
             if n > sample else np.arange(n))
    rv = np.ones(len(cells))
    for k, i in enumerate(cells):
        r, _ = stats.spearmanr(mygetidx(matrix, i, sparse=sparse, axis=0),
                               mygetidx(raw, i, axis=0))
        rv[k] = r
    return rv

def plot_mono(matrix, raw, ax, mono_vals=None):
    # mono_vals may be precomputed (PFlog is rank-monotone by construction, so
    # every per-cell Spearman is 1.0 -- no densify / Spearman loop needed).
    x = mono(matrix, raw) if mono_vals is None else mono_vals
    p = {"xlabel": "Spearman r", "ylabel": "Frequency", "xlim": (-1.2, 1.2)}
    close = np.all(np.allclose(x, x[0]))
    if close:
        weights = np.ones(len(x)) / len(x)
        x = np.array([1] * len(x))
        ax.hist(x, facecolor=colors["mono"], edgecolor="k", weights=weights)
    else:
        weights = np.ones(len(x)) / len(x)
        ax.hist(x, facecolor=colors["mono"], weights=weights, edgecolor="k")
    xmean = x.mean()
    ax.axvline(xmean, linestyle="--", color="darkgray", label=f"mean: {xmean:,.2f}")
    ax.set(**p)
    ax.legend(prop={"size": 12})
    return (ax, xmean)

def clr_gene_var(sclr):
    # Per-gene variance of the dense shifted-CLR  Z = L - m.1^T  computed from the
    # sparse log1p(PF) matrix L and the per-cell center vector m WITHOUT densifying:
    #   Var_g(Z) = Var_g(L) + Var(m) - 2 Cov_g(L, m).
    L = sclr.sparse
    m = np.asarray(sclr.row_center, dtype=float).ravel()
    n = L.shape[0]
    EL = np.asarray(L.sum(0)).ravel() / n
    VarL = np.asarray(L.multiply(L).sum(0)).ravel() / n - EL**2
    Cov = np.asarray(L.T.dot(m)).ravel() / n - EL * m.mean()
    return VarL + m.var() - 2 * Cov

def read_data(base_data_fn, max_cells=25000, seed=0):
    from metrics_matrix import compute_overdispersion
    from scclr_pflog import normalize_pflog
    data = {}

    for title in mtx_labels:
        in_fn = os.path.join(base_data_fn, f"{title}.mtx.gz")
        data[txlabel[title]] = mmread(in_fn).tocsr()

    # Dense methods loaded from csv.gz. Missing files are silently skipped so
    # plots still produce on datasets that haven't computed all methods yet.
    for title in ["sctransform", "cp10k_log_scale"]:
        in_fn = os.path.join(base_data_fn, f"{title}.csv.gz")
        if os.path.exists(in_fn):
            data[txlabel[title]] = pd.read_csv(in_fn, header=None, compression="gzip").values

    # Cap cells (the SAME subset across every method) on very large datasets. PFlog
    # is now sparse (scclr, below), but sctransform and scalelog1pCP10k are inherently
    # dense (loaded from multi-GB csv.gz), so the cap still bounds their memory under
    # parallel batch rendering. A seeded subsample is statistically indistinguishable
    # for these per-gene / per-cell summaries.
    raw = data[txlabel["raw"]]
    if raw.shape[0] > max_cells:
        idx = np.sort(np.random.default_rng(seed).choice(raw.shape[0], max_cells, replace=False))
        for k, v in list(data.items()):
            data[k] = v[idx]
        raw = data[txlabel["raw"]]

    # PFlog (shift. CLR) is current scclr PFlog:
    # center(log1p(4*alpha*x)).  scclr stores it as sparse log terms plus a
    # per-cell center vector, so the panel is memory-safe even on the largest
    # matrices. plot_data derives panel stats from the sparse-plus-rank-one form.
    alpha = float(compute_overdispersion(raw))
    data["PFlog (shift. CLR)"] = normalize_pflog(raw, alpha=alpha)
    return data

# 'pf_log_pf' intentionally NOT in mtx_labels: PFlog is computed on the fly as
# current scclr PFlog (see read_data), not read from the legacy
# pf_log_pf.mtx.gz artifact.
mtx_labels = ['raw', 'pf', 'log', 'pf_log', 'cpm_log', 'cp10k_log', "sqrt"]

# "PFlog (shift. CLR)" = current scclr shifted CLR, computed in read_data.
labels = [
    'raw',
     'PF',
     "sqrt",
     'log1p',
     'log1pCP10k',
     'log1pCPM',
     'scalelog1pCP10k',
     'sctransform',
     'log1pPF',
     'PFlog (shift. CLR)',
]

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
  'pf_log_pf': 'legacy PF-log-PF',
}

def setup_plot(ds, shape):
    # 3 subplots per method (variance stabilization + depth normalization + monotonicity).
    n_methods = len(labels)
    n_cols = 2
    n_rows = (n_methods + n_cols - 1) // n_cols
    fig = plt.figure(figsize=(6*3, n_rows*3))
    fig.suptitle(fr"{ds} ({shape[0]:,.0f} $\times$ {shape[1]:,.0f})", y=0.92)

    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig, wspace=0.15, hspace=0.75)
    axs = []
    for i in range(n_rows):
        for j in range(n_cols):
            ig = gs[i,j].subgridspec(1, 3, wspace=0.4)
            ax1 = fig.add_subplot(ig[0, 0])
            ax2 = fig.add_subplot(ig[0, 1])
            ax3 = fig.add_subplot(ig[0, 2])
            axs.append((ax1, ax2, ax3))
    return (fig, axs)

_CVQ_CACHE = None
def _cvq_lookup():
    """Per-dataset, per-method CV(q) from the committed Fig 1b CSVs, used for the
    variance-stabilization panel annotation (functional transforms = model-based
    delta-method q; sctransform = empirical CV of full-vst residual gene variances)."""
    global _CVQ_CACHE
    if _CVQ_CACHE is not None:
        return _CVQ_CACHE
    import pandas as pd
    figdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis", "figures")
    look = {}
    try:
        f = pd.read_csv(os.path.join(figdir, "cvq_per_method_subset_filtered.csv"))
        cols = ["raw","PF","sqrt","log1p","log1pCP10k","log1pCPM","scalelog1pCP10k","log1pPF","PFlog"]
        for _, r in f.iterrows():
            d = {c: float(r[c]) for c in cols if c in r and r[c] == r[c]}
            if "PFlog" in d:
                d["PFlog (shift. CLR)"] = d.pop("PFlog")
            look[r["ds"]] = d
    except Exception as e:
        print(f"[plot_all] CV(q) functional CSV not loaded: {e}", file=sys.stderr)
    try:
        s = pd.read_csv(os.path.join(figdir, "sct_cvq_filtered.csv"))
        for _, r in s.iterrows():
            if r["cv_filt"] == r["cv_filt"]:
                look.setdefault(r["ds"], {})["sctransform"] = float(r["cv_filt"])
    except Exception as e:
        print(f"[plot_all] CV(q) sctransform CSV not loaded: {e}", file=sys.stderr)
    _CVQ_CACHE = look
    return look

def plot_data(axs, data, ds=""):
    realds = ds[:-len("_subset_genes")] if ds.endswith("_subset_genes") else ds
    cvq = _cvq_lookup().get(realds, {})
    raw = data["raw"]
    raw_genevar = myvar(raw, axis=0)
    raw_genemean = mymean(raw, axis=0)
    raw_cellsum  = mysum(raw, axis=1)

    minlim = min(np.min(raw_genevar), np.min(raw_genemean)) * 0.1
    maxlim = max(np.max(raw_genevar), np.max(raw_genemean)) * 10
    metrics = defaultdict(dict)
    for (ax1, ax2, ax3), title in zip(axs, labels):
        m = data.get(title)
        if m is None:
            ax1.set_visible(False)
            ax2.axis("off")
            ax3.set_visible(False)
            ax2.text(0.5, 0.5, 'not computed',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax2.transAxes)
        else:
            try:
                cq = cvq.get(title)   # CV(q) for this method (None -> falls back to cov_gene)
                if title == "PFlog (shift. CLR)":
                    # scclr sparse shifted-CLR (ShiftedCLR object): derive panel stats
                    # from the sparse log1p(PF) + per-cell center -- never densified.
                    # Centered => per-cell sums are 0 (depth removed); rank-monotone by
                    # construction => every per-cell Spearman is 1.
                    n = m.shape[0]
                    (_, cov_gene) = plot_meanvar(None, raw_genemean, minlim=minlim, maxlim=maxlim,
                                                 ax=ax1, gvar=clr_gene_var(m), cq=cq)
                    (_, r2_depth) = plot_depth(None, raw_cellsum, ax2, cell_sums=np.zeros(n))
                    (_, r_mono)   = plot_mono(None, raw, ax3, mono_vals=np.ones(n))
                else:
                    (_, cov_gene) = plot_meanvar(m, raw_genemean, minlim = minlim, maxlim = maxlim, ax=ax1, cq=cq)
                    (_, r2_depth) = plot_depth(m, raw_cellsum, ax2)
                    (_, r_mono) = plot_mono(m, raw, ax3)
            except (MemoryError, KeyboardInterrupt):
                # never silently publish a blank panel on a resource failure
                raise
            except Exception as e:
                print(f"[plot_all] {ds}: '{title}' panel failed: {type(e).__name__}: {e}",
                      file=sys.stderr)
                ax1.set_visible(False)
                ax2.axis("off")
                ax3.set_visible(False)
                ax2.text(0.5, 0.5, 'NaN',
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax2.transAxes)

        disp = "sctransform v2" if title == "sctransform" else title
        ax2.set_title(disp, fontsize=20, weight="bold")

    return

def main(ds, base_data_fn, base_out_fn):
#     ds = "GSM3711776"
#     base_data_fn =  os.path.join("synchromesh/data/", ds)
#     base_out_fn =  os.path.join("synchromesh/data/", ds)
    data = read_data(base_data_fn)
    # mmread returns COO; per-cell getrow() in mono() needs CSR to avoid
    # COO->LIL conversion on every call. Single conversion here, ~10x speedup on big sets.
    for k, v in data.items():
        if issparse(v):
            data[k] = v.tocsr()

    shape = data["raw"].shape
    fig, axs = setup_plot(ds, shape)
    plot_data(axs, data, ds=ds)

    fig.savefig(os.path.join(base_out_fn, f"{ds}_method_comparison.pdf"), facecolor='white', transparent=False, dpi=300, bbox_inches="tight")
    return

if __name__ == "__main__":
    # ds, base_data_fn, base_out_fn
    main(sys.argv[1], sys.argv[2], sys.argv[3])
