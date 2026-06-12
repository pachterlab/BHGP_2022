from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd


matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
FIGURE_OUT = OUT.parents[1] / "output"


def main():
    df_long = pd.read_csv(OUT / "smartseq3_fibroblasts_procrustes_10d_displacements.csv")
    logpf = df_long[df_long["method"] == "log1pPF"]["displacement"].to_numpy()
    pflogpf = df_long[df_long["method"] == "PFlog1pPF"]["displacement"].to_numpy()
    if len(logpf) != len(pflogpf):
        raise ValueError(f"unequal lengths: {len(logpf)} log1pPF, {len(pflogpf)} PFlog1pPF")

    paired = pd.DataFrame(
        {
            "cell_index": np.arange(len(logpf)),
            "log1pPF_distance": logpf,
            "PFlog1pPF_distance": pflogpf,
            "ratio_PFlogPF_over_logPF": pflogpf / np.maximum(logpf, 1e-15),
        }
    )
    paired.to_csv(OUT / "smartseq3_fibroblasts_paired_distance_scatter.csv", index=False)

    frac_better = float(np.mean(pflogpf < logpf))
    med_logpf = float(np.median(logpf))
    med_pflogpf = float(np.median(pflogpf))
    mean_logpf = float(np.mean(logpf))
    mean_pflogpf = float(np.mean(pflogpf))
    lim = float(max(np.quantile(logpf, 0.995), np.quantile(pflogpf, 0.995)) * 1.08)

    fig, ax = plt.subplots(figsize=(3.45, 3.2), constrained_layout=True)
    ax.scatter(logpf, pflogpf, s=12, color="#2B6CB0", alpha=0.42, linewidths=0)
    ax.plot([0, lim], [0, lim], color="#444444", lw=1.0, ls="--", zorder=0)
    ax.axvline(med_logpf, color="#2B6CB0", lw=1.1, ls=":")
    ax.axhline(med_pflogpf, color="#D55E00", lw=1.1, ls=":")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("log1pPF paired distance")
    ax.set_ylabel("PFlog1pPF paired distance")
    ax.text(
        0.04,
        0.96,
        f"{frac_better:.1%} below diagonal\nmedian {med_logpf:.4f} vs {med_pflogpf:.4f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    FIGURE_OUT.mkdir(exist_ok=True)
    fig.savefig(FIGURE_OUT / "smartseq3_fibroblasts_paired_distance_scatter.png", dpi=300)
    fig.savefig(FIGURE_OUT / "smartseq3_fibroblasts_paired_distance_scatter.pdf")

    summary = pd.DataFrame(
        [
            {
                "n_cells": len(logpf),
                "fraction_PFlogPF_less_than_logPF": frac_better,
                "median_log1pPF_distance": med_logpf,
                "median_PFlog1pPF_distance": med_pflogpf,
                "mean_log1pPF_distance": mean_logpf,
                "mean_PFlog1pPF_distance": mean_pflogpf,
            }
        ]
    )
    summary.to_csv(OUT / "smartseq3_fibroblasts_paired_distance_scatter_summary.csv", index=False)
    print(summary.to_string(index=False))
    print("outputs", OUT)


if __name__ == "__main__":
    main()
