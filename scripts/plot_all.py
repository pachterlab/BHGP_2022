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

cividis = matplotlib.cm.get_cmap("cividis")
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
        return (a_squared.mean(axis) - np.square(a.mean(axis))).A.ravel()
    else:
        return np.var(a, axis=axis)
    
def mymean(a, sparse=True, axis=None):
    if sparse:
        return a.mean(axis=axis).A.ravel()
    else:
        return a.mean(axis=axis)
    
def mysum(a, sparse=True, axis=None):
    if sparse:
        return a.sum(axis=axis).A.ravel()
    else:
        return a.sum(axis=axis)
    
def mygetidx(a, idx, sparse=True, axis=0):
    if sparse:
        if axis == 0:
            return a.getrow(idx).A.ravel()
        elif axis == 1:
            return a.getcol(idx).A.ravel()
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

def mono(matrix, raw):
    sparse = issparse(matrix)
    rv = np.ones(matrix.shape[0])
    if sparse:
        # sparse matrices are all monotonic, you can run this to verify, will return
        # 1's for all cell-cell comparisons (spearman)
        return rv
    for i in range(matrix.shape[0]):
        r, p  = stats.spearmanr(mygetidx(matrix, i, sparse=sparse, axis=0), mygetidx(raw, i, axis=0))
        rv[i] = r
    return rv

def plot_meanvar(mtx, raw_mean, minlim = 1e-4, maxlim = 1e5, ax=None):
    p = {
        "xlabel": "Gene mean",
        "ylabel": "Gene variance",
        "xscale": "log",
        "yscale": "log",
        "xlim": (minlim, maxlim),
    }
    
    sparse = issparse(mtx)
    gvar = myvar(mtx, sparse=sparse, axis=0)
    gcov = np.sqrt(np.var(gvar))/np.mean(gvar)

    y = gvar
    yy = (y-y.mean())/np.sqrt(np.var(y))

    ax.scatter(raw_mean, y, facecolor=colors["gene"], alpha=alpha, edgecolor="k", label=f"CoV: {gcov:,.1f}")
    ax.legend(prop={"size": 12})
    ax.set(**p)
    yex(ax)
    return (ax, gcov)

def plot_depth(mtx, raw_cell_counts, ax):
    x = raw_cell_counts
    sparse = issparse(mtx)
    y = mysum(mtx, sparse=sparse, axis=1)
    
    minx, maxx = min(x), max(x)
    miny, maxy = min(y), max(y)
    maxy = maxy - miny

    xx = (x - minx)/maxx
    yy = (y - miny)/maxy
    
    close = np.all(np.allclose(y, y[0]))
    if close:
        yy = [1]*len(y)
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

def plot_mono(matrix, raw, ax):
    x = mono(matrix, raw)
    p = {
        "xlabel": "Spearman r",
        "ylabel": "Frequency",
        "xlim": (-1.2, 1.2)
    }
    close = np.all(np.allclose(x, x[0]))
    if close:
        weights=np.ones(len(x)) / len(x)
        x = np.array([1] * len(x))
        ax.hist(x, facecolor=colors["mono"], edgecolor="k", weights=weights)
    else:
        weights=np.ones(len(x)) / len(x)
        ax.hist(x, facecolor=colors["mono"], weights=weights, edgecolor="k")
    xmean = x.mean()
    ax.axvline(x.mean(), linestyle="--", color="darkgray", label=f"mean: {xmean:,.2f}")
    ax.set(**p)
    ax.legend(prop={"size": 12})
    return (ax, xmean)

def read_data(base_data_fn):
    data = {}

    for title in mtx_labels:
        in_fn = os.path.join(base_data_fn, f"{title}.mtx.gz")
        data[txlabel[title]] = mmread(in_fn)

    title = "sctransform"
    in_fn = os.path.join(base_data_fn, f"{title}.csv.gz")
    data[txlabel[title]] = pd.read_csv(in_fn, header=None, compression="gzip").values

    title = "cp10k_log_scale"
    in_fn = os.path.join(base_data_fn, f"{title}.csv.gz")
    data[txlabel[title]] = pd.read_csv(in_fn, header=None, compression="gzip").values
    return data

mtx_labels = ['raw', 'pf', 'log', 'pf_log', 'pf_log_pf', 'cpm_log', 'cp10k_log', "sqrt"]

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
     'PFlog1pPF'
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
  'pf_log_pf': 'PFlog1pPF'
}

def setup_plot(ds, shape):
    fig = plt.figure(figsize=(6*3,5*3))
    fig.suptitle(fr"{ds} ({shape[0]:,.0f} $\times$ {shape[1]:,.0f})", y=0.92)

    gs = gridspec.GridSpec(5, 2, figure=fig, wspace=0.15, hspace=0.75)
    axs = []
    for i in range(5):
        for j in range(2):
            ig = gs[i,j].subgridspec(1, 3, wspace=0.4)
            ax1 = fig.add_subplot(ig[0, 0])
            ax2 = fig.add_subplot(ig[0, 1])
            ax3 = fig.add_subplot(ig[0, 2])
            axs.append((ax1, ax2, ax3))
    return (fig, axs)

def plot_data(axs, data):
    raw = data["raw"]
    raw_genevar = myvar(raw, axis=0)
    raw_genemean = mymean(raw, axis=0)
    raw_cellsum  = mysum(raw, axis=1)
    
    minlim = min(np.min(raw_genevar), np.min(raw_genemean)) * 0.1
    maxlim = max(np.max(raw_genevar), np.max(raw_genemean)) * 10
    metrics = defaultdict(dict)
    for (ax1, ax2, ax3), title in zip(axs, labels):
        m = data[title]
        try:
            (_, cov_gene) = plot_meanvar(m, raw_genemean, minlim = minlim, maxlim = maxlim, ax=ax1)
            (_, r2_depth) = plot_depth(m, raw_cellsum, ax2)
            (_, r_mono) = plot_mono(m, raw, ax3)
        except:
            ax1.set_visible(False)
            ax2.axis("off")
            ax3.set_visible(False)
            ax2.text(0.5, 0.5, 'NaN',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax2.transAxes)

        ax2.set_title(title, fontsize=20, weight="bold")
    
    return

def main(ds, base_data_fn, base_out_fn):
#     ds = "GSM3711776"
#     base_data_fn =  os.path.join("synchromesh/data/", ds)
#     base_out_fn =  os.path.join("synchromesh/data/", ds)
    data = read_data(base_data_fn)
    
    shape = data["raw"].shape
    fig, axs = setup_plot(ds, shape)
    plot_data(axs, data)

    fig.savefig(os.path.join(base_out_fn, f"{ds}_method_comparison.pdf"), facecolor='white', transparent=False, dpi=300, bbox_inches="tight")
    return

if __name__ == "__main__":
    # ds, base_data_fn, base_out_fn
    main(sys.argv[1], sys.argv[2], sys.argv[3])
