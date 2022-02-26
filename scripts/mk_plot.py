#!/usr/bin/env python3

from operator import index
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.io import mmread
import pandas as pd

matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 15
import numpy as np

alpha = 0.25


def yex(ax):
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, "k-", alpha=0.75, zorder=0)
    ax.set_aspect("equal")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    return ax


def plot_meanvar(mtx, raw_mean, minlim=1e-4, maxlim=1e5, ax=None):
    p = {
        "xlabel": "Gene mean",
        "ylabel": "Gene variance",
        "xscale": "log",
        "yscale": "log",
        "xlim": (minlim, maxlim),
    }

    y = np.var(mtx, axis=0)

    ax.scatter(raw_mean, y, facecolor="#D43F3A", alpha=alpha, edgecolor="k")
    ax.axhline(1, color="k", linestyle="--")
    ax.set(**p)
    yex(ax)
    return ax


def plot_depth(mtx, raw_cell_counts, ax):
    x = np.sum(mtx, axis=1)

    mean = x.mean()
    stdev = np.sqrt(np.var(raw_cell_counts))

    close = np.all(np.allclose(x, x[0]))

    if close:
        weights = np.ones(len(x)) / len(x)
        ax.hist([0] * len(x), facecolor="#D43F3A", weights=weights, edgecolor="k")
    else:
        weights = np.ones(len(x)) / len(x)
        xx = (x - mean) / stdev
        ax.hist(xx, facecolor="#D43F3A", weights=weights, edgecolor="k")
    p = {"xlabel": "Z-score (cell counts)", "ylabel": "Frequency", "xlim": (-3, 3)}
    ax.set(**p)
    return ax


def mono(matrix, raw):
    rv = np.zeros(matrix.shape[0])
    for i in range(matrix.shape[0]):
        r, p = stats.spearmanr(matrix[i], raw[i])
        rv[i] = r
    return rv


def plot_mono(matrix, raw, ax):
    x = mono(matrix, raw)
    p = {"xlabel": "Spearman r", "ylabel": "Frequency", "xlim": (-0.1, 1.2)}
    close = np.all(np.allclose(x, x[0]))
    if close:
        weights = np.ones(len(x)) / len(x)
        ax.hist([1] * len(x), facecolor="#D43F3A", edgecolor="k", weights=weights)
    else:
        weights = np.ones(len(x)) / len(x)
        ax.hist(x, facecolor="#D43F3A", weights=weights, edgecolor="k")
    ax.set(**p)
    return ax


def read_data(base_data_fn):
    data = {}
    titles = ["raw", "pf", "log", "pf_log", "pf_log_pf", "cpm_log", "cp10k_log"]
    for title in titles:
        print(f"loading {title}")
        in_fn = os.path.join(base_data_fn, f"{title}.mtx")
        data[title] = mmread(in_fn).toarray()

    title = "sctransform"
    print(f"loading {title}")
    in_fn = os.path.join("./", f"{title}.csv")
    data[title] = pd.read_csv(in_fn, header=None).values

    title = "cp10k_log_scale"
    print(f"loading {title}")
    in_fn = os.path.join("./", f"{title}.csv")
    data[title] = pd.read_csv(in_fn, header=None).values
    return data


def main(ds, base_data_fn, base_out_fn):
    data = read_data(base_data_fn)

    print(f"saving {ds}.png")
    fig = plt.figure(figsize=(6 * 3, 5 * 3))
    fig.suptitle(ds, y=0.92)
    gs = gridspec.GridSpec(5, 2, figure=fig, wspace=0.15, hspace=0.75)
    axs = []
    for i in range(5):
        for j in range(2):
            ig = gs[i, j].subgridspec(1, 3, wspace=0.4)
            ax1 = fig.add_subplot(ig[0, 0])
            ax2 = fig.add_subplot(ig[0, 1])
            ax3 = fig.add_subplot(ig[0, 2])
            axs.append((ax1, ax2, ax3))

    labels = [
        "raw",
        "pf",
        "log",
        "pf_log",
        "pf_log_pf",
        "cpm_log",
        "cp10k_log",
        "cp10k_log_scale",
        "sctransform",
    ]
    raw = data["raw"]
    ax1, ax2, ax3 = axs.pop()
    ax1.remove()
    ax2.remove()
    ax3.remove()

    minlim = min(np.min(np.var(raw, 0)), np.min(np.mean(raw, 0))) * 0.1
    maxlim = max(np.max(np.var(raw, 0)), np.max(np.mean(raw, 0))) * 10
    for (ax1, ax2, ax3), title in zip(axs, labels):
        m = data[title]
        plot_meanvar(m, raw.mean(0), minlim=minlim, maxlim=maxlim, ax=ax1)
        plot_depth(m, raw.sum(1), ax2)
        plot_mono(m, raw, ax3)
        ax2.set_title(title, fontsize=20, weight="bold")

    out_fn = os.path.join(base_out_fn, f"{ds}.png")
    fig.savefig(
        out_fn, facecolor="white", transparent=False, dpi=200, bbox_inches="tight"
    )
    fig.show()


if __name__ == "__main__":
    # ds, base_data_fn, base_out_fn
    main(sys.argv[1], sys.argv[2], sys.argv[3])
