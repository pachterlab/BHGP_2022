#!/usr/bin/env python3
"""Plot depth confounding diagnostics for the linear/random-walk simulations.

The Ahlmann-Eltze--Huber linear-walk and random-walk simulators use real
cell-specific depths when sampling counts. This diagnostic asks whether those
depths are structured with respect to the simulator ground-truth neighborhoods.
If ground-truth neighbors are much closer in depth than random cell pairs, then
technical depth is predictive of the benchmark target.
"""

from __future__ import annotations

from pathlib import Path
import importlib.util

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = (ROOT.parent / "figures").resolve()
CACHE = FIGURE_DIR / "scclr_cache"
REPO_OUTPUT = ROOT / "PFlogPF" / "output"
OUT_PNG = REPO_OUTPUT / "walk_depth_confounding.png"
OUT_PDF = REPO_OUTPUT / "walk_depth_confounding.pdf"
OUT_TSV = REPO_OUTPUT / "walk_depth_confounding_summary.tsv"

SIMULATORS = {
    "linear_walk": "Linear walk",
    "random_walk": "Random walk",
}


def ensure_walk_cache(simulator: str, seed: int) -> None:
    counts = CACHE / f"simulation_allv2_{simulator}_seed{seed}_counts.mtx"
    truth_knn = CACHE / f"simulation_allv2_{simulator}_seed{seed}_truth_knn.tsv"
    if counts.exists() and truth_knn.exists():
        return

    helper_path = ROOT / "scripts" / "recompute_fig2abc_scclr.py"
    spec = importlib.util.spec_from_file_location("recompute_fig2abc_scclr", helper_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {helper_path}")
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    helpers.ensure_dirs()
    helpers.generate_simulation_mtx(simulator, seed)


def read_counts(simulator: str, seed: int):
    ensure_walk_cache(simulator, seed)
    path = CACHE / f"simulation_allv2_{simulator}_seed{seed}_counts.mtx"
    if not path.exists():
        raise FileNotFoundError(path)
    return mmread(path).tocsr()


def read_truth_knn(simulator: str, seed: int) -> np.ndarray:
    ensure_walk_cache(simulator, seed)
    path = CACHE / f"simulation_allv2_{simulator}_seed{seed}_truth_knn.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    knn = np.loadtxt(path, delimiter="\t", dtype=np.int64)
    knn = np.atleast_2d(knn)
    if knn.min() == 1:
        knn = knn - 1
    return knn


def rolling_median(y: np.ndarray, window: int = 151) -> tuple[np.ndarray, np.ndarray]:
    half = window // 2
    x_out = []
    y_out = []
    for center in range(half, len(y) - half, max(1, window // 8)):
        x_out.append(center)
        y_out.append(np.median(y[(center - half) : (center + half + 1)]))
    return np.asarray(x_out), np.asarray(y_out)


def random_knn(n_cells: int, k: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty((n_cells, k), dtype=np.int64)
    universe = np.arange(n_cells)
    for i in range(n_cells):
        pool = universe[universe != i]
        out[i] = rng.choice(pool, size=k, replace=False)
    return out


def mean_depth_gap(depth: np.ndarray, knn: np.ndarray) -> float:
    rows = np.arange(knn.shape[0])[:, None]
    return float(np.mean(np.abs(depth[rows] - depth[knn])))


def collect_summary() -> pd.DataFrame:
    rows = []
    for simulator in SIMULATORS:
        for seed in range(1, 6):
            counts = read_counts(simulator, seed)
            depth = np.asarray(counts.sum(axis=0)).ravel()
            truth = read_truth_knn(simulator, seed)
            random = random_knn(len(depth), truth.shape[1], seed=1000 + seed)
            rho, pvalue = spearmanr(np.arange(len(depth)), depth)
            rows.append(
                {
                    "simulator": simulator,
                    "seed": seed,
                    "n_cells": len(depth),
                    "median_depth": float(np.median(depth)),
                    "cell_order_spearman_r": float(rho),
                    "cell_order_spearman_p": float(pvalue),
                    "truth_neighbor_depth_gap": mean_depth_gap(depth, truth),
                    "random_pair_depth_gap": mean_depth_gap(depth, random),
                }
            )
    out = pd.DataFrame(rows)
    out["truth_to_random_gap_ratio"] = (
        out["truth_neighbor_depth_gap"] / out["random_pair_depth_gap"]
    )
    return out


def plot(summary: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
        }
    )

    colors = {
        "linear_walk": "#2B6CB0",
        "random_walk": "#C05621",
        "Truth neighbors": "#2F855A",
        "Random pairs": "#8A8F98",
    }

    fig = plt.figure(figsize=(7.1, 5.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0])
    ax_linear = fig.add_subplot(gs[0, 0])
    b_gs = gs[0, 1].subgridspec(2, 1, height_ratios=[0.36, 0.64], hspace=0.05)
    ax_random_top = fig.add_subplot(b_gs[0, 0])
    ax_random = fig.add_subplot(b_gs[1, 0], sharex=ax_random_top)
    ax = fig.add_subplot(gs[1, :])

    counts = read_counts("linear_walk", 1)
    depth = np.asarray(counts.sum(axis=0)).ravel()
    x = np.arange(1, len(depth) + 1)
    ax_linear.scatter(x, depth, s=5, alpha=0.22, color=colors["linear_walk"], linewidths=0)
    rx, ry = rolling_median(depth)
    ax_linear.plot(rx + 1, ry, color="black", lw=1.4)
    ax_linear.set_title("Linear walk", loc="left", fontweight="bold")
    ax_linear.set_xlabel("Simulator cell order")
    ax_linear.set_ylabel("Total counts")
    ax_linear.text(
        -0.12,
        1.06,
        "a",
        transform=ax_linear.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
    )

    counts = read_counts("random_walk", 1)
    depth = np.asarray(counts.sum(axis=0)).ravel()
    x = np.arange(1, len(depth) + 1)
    for target_ax in (ax_random_top, ax_random):
        target_ax.scatter(x, depth, s=5, alpha=0.22, color=colors["random_walk"], linewidths=0)
    rx, ry = rolling_median(depth)
    ax_random.plot(rx + 1, ry, color="black", lw=1.4)
    ax_random_top.set_ylim(70000, 620000)
    ax_random.set_ylim(0, 3500)
    ax_random_top.set_title("Random walk", loc="left", fontweight="bold")
    ax_random.set_xlabel("Simulator cell order")
    ax_random.set_ylabel("Total counts")
    ax_random_top.spines["bottom"].set_visible(False)
    ax_random.spines["top"].set_visible(False)
    ax_random_top.tick_params(labelbottom=False, bottom=False)
    ax_random.tick_params(top=False)
    ax_random_top.text(
        -0.12,
        1.10,
        "b",
        transform=ax_random_top.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
    )
    d = 0.012
    kwargs = dict(color="black", clip_on=False, lw=0.8)
    ax_random_top.plot((-d, +d), (-d, +d), transform=ax_random_top.transAxes, **kwargs)
    ax_random_top.plot((1 - d, 1 + d), (-d, +d), transform=ax_random_top.transAxes, **kwargs)
    ax_random.plot((-d, +d), (1 - d, 1 + d), transform=ax_random.transAxes, **kwargs)
    ax_random.plot((1 - d, 1 + d), (1 - d, 1 + d), transform=ax_random.transAxes, **kwargs)

    ax.text(
        -0.035,
        1.07,
        "c",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
    )
    x_positions = []
    labels = []
    width = 0.32
    for i, simulator in enumerate(SIMULATORS):
        sub = summary[summary["simulator"] == simulator]
        x0 = i * 1.35
        truth_vals = sub["truth_neighbor_depth_gap"].to_numpy()
        rand_vals = sub["random_pair_depth_gap"].to_numpy()
        for offset, vals, label, color in [
            (-width / 2, truth_vals, "Truth neighbors", colors["Truth neighbors"]),
            (width / 2, rand_vals, "Random pairs", colors["Random pairs"]),
        ]:
            xpos = x0 + offset
            ax.bar(
                xpos,
                vals.mean(),
                width=width,
                color=color,
                edgecolor="black",
                linewidth=0.5,
                label=label if i == 0 else None,
            )
            ax.errorbar(
                xpos,
                vals.mean(),
                yerr=vals.std(ddof=1),
                color="black",
                capsize=2,
                lw=0.8,
            )
            jitter = np.linspace(-0.045, 0.045, len(vals))
            ax.scatter(
                np.full(len(vals), xpos) + jitter,
                vals,
                s=18,
                color="white",
                edgecolor="black",
                linewidth=0.5,
                zorder=3,
            )
        ratio = sub["truth_to_random_gap_ratio"].mean()
        ax.text(
            x0,
            max(truth_vals.mean(), rand_vals.mean()) + 0.035,
            f"ratio = {ratio:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
        x_positions.append(x0)
        labels.append(SIMULATORS[simulator])

    ax.set_xticks(x_positions, labels)
    ax.set_ylabel("Mean absolute depth difference")
    ax.set_title("Ground-truth neighbors are depth-similar", loc="left", fontweight="bold")
    ax.legend(
        frameon=False,
        ncol=2,
        loc="upper right",
        bbox_to_anchor=(1.03, 1.17),
        fontsize=7.5,
        handlelength=1.2,
        columnspacing=1.1,
    )
    ax.set_ylim(bottom=0)

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")


def main() -> None:
    REPO_OUTPUT.mkdir(parents=True, exist_ok=True)
    summary = collect_summary()
    summary.to_csv(OUT_TSV, sep="\t", index=False)
    plot(summary)
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_PDF}")
    print(f"wrote {OUT_TSV}")
    print(summary.groupby("simulator")["truth_to_random_gap_ratio"].agg(["mean", "std"]))


if __name__ == "__main__":
    main()
