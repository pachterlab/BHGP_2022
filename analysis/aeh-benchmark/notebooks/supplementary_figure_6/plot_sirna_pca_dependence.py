#!/usr/bin/env python3

from pathlib import Path
import os
import tempfile

tmp_dir = Path(tempfile.gettempdir())
os.environ.setdefault("MPLCONFIGDIR", str(tmp_dir / "mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(tmp_dir))

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PF_ROOT = SCRIPT_DIR.parents[1]
OUT_DIR = PF_ROOT / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RESULTS = (
    PF_ROOT
    / "benchmark"
    / "output"
    / "benchmark_results"
    / "downsampling_results.tsv"
)
CLR_OVERRIDE = SCRIPT_DIR / "sirna_pca_dependence_clr_alpha_k.tsv"

FAMILY_MAP = {
    "logp1": "delta_method",
    "acosh": "delta_method",
    "logp_alpha": "delta_method",
    "logp_cpm": "delta_method",
    "logp1_size_normed": "delta_method",
    "logp1_hvg": "delta_method",
    "logp1_zscore": "delta_method",
    "logp1_hvg_zscore": "delta_method",
    "clr": "compositional",
    "clr_alpha": "compositional",
    "pearson_clip": "glm_residual",
    "sctransform": "glm_residual",
    "pearson_analytic": "glm_residual",
    "rand_quantile": "glm_residual",
    "pearson": "glm_residual",
    "pearson_clip_hvg": "glm_residual",
    "pearson_clip_zscore": "glm_residual",
    "pearson_clip_hvg_zscore": "glm_residual",
    "sanity_map": "latent_expr",
    "sanity_dists": "latent_expr",
    "dino": "latent_expr",
    "normalisr_norm": "latent_expr",
    "scgpt": "value_binning",
    "glmpca": "count_model",
    "newwave": "count_model",
    "raw_counts": "negative_control",
    "scaled_raw_counts": "negative_control",
}

FAMILY_COLORS = {
    "delta_method": "#66C2A5",
    "compositional": "#E41A1C",
    "glm_residual": "#FC8D62",
    "latent_expr": "#8DA0CB",
    "value_binning": "#A6D854",
    "count_model": "#e78ac3",
    "negative_control": "#7d7d7d",
}

DIRECT_LABELS = {
    "logp1": "log1pPF",
    "clr": "PFlog",
}


def main() -> None:
    results = pd.read_csv(RESULTS, sep="\t")
    if CLR_OVERRIDE.exists():
        overrides = pd.read_csv(CLR_OVERRIDE, sep="\t")
        key_cols = ["dataset", "seed", "pca_dim", "knn"]
        override_keys = overrides[key_cols].drop_duplicates()
        stale = results.merge(override_keys.assign(_replace=True), how="left", on=key_cols)
        results = pd.concat(
            [
                stale[~((stale["transformation"] == "clr") & stale["_replace"].fillna(False))]
                .drop(columns=["_replace"]),
                overrides,
            ],
            ignore_index=True,
            sort=False,
        )
    else:
        raise FileNotFoundError(
            f"Missing {CLR_OVERRIDE}. Run recompute_sirna_pca_dependence_scclr.py first."
        )

    sirna = results[
        (results["dataset"] == "smartSeq3_siRNA_knockdown")
        & (results["knn"] == 50)
        & (results["pca_dim"].isin([5, 10, 50, 100]))
        & (results["alpha"].astype(str).isin(["True", "False", "TRUE", "FALSE"]))
    ]

    plot_data = (
        sirna.groupby(["transformation", "pca_dim"], as_index=False)["overlap"]
        .mean()
        .rename(columns={"overlap": "knn_overlap"})
        .sort_values(["transformation", "pca_dim"])
    )

    source_data = SCRIPT_DIR / "sirna_pca_dependence_source_data.csv"
    plot_data.to_csv(source_data, index=False)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
        }
    )

    fig, ax = plt.subplots(figsize=(110 / 25.4, 75 / 25.4))
    ax.set_xscale("log")
    ax.set_xlim(4.6, 145)
    ax.set_ylim(0, 37)
    ax.set_xticks([5, 10, 50, 100])
    ax.set_xticklabels(["5", "10", "50", "100"])
    ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35])
    ax.grid(axis="y", color="#e6e6e6", linewidth=0.8)
    ax.grid(axis="x", color="#e6e6e6", linewidth=0.8)

    for trans, trans_data in plot_data.groupby("transformation"):
        family = FAMILY_MAP.get(trans, "negative_control")
        color = FAMILY_COLORS[family]
        ax.plot(
            trans_data["pca_dim"],
            trans_data["knn_overlap"],
            color=color,
            linewidth=1.5 if trans == "clr" else 0.8,
            linestyle="--" if trans == "sanity_dists" else "-",
            solid_capstyle="round",
            zorder=3 if trans == "clr" else 2,
        )

    for trans, label in DIRECT_LABELS.items():
        endpoint = plot_data[
            (plot_data["transformation"] == trans) & (plot_data["pca_dim"] == 100)
        ].iloc[0]
        family = FAMILY_MAP[trans]
        ax.text(
            106,
            endpoint["knn_overlap"],
            label,
            color=FAMILY_COLORS[family],
            va="center",
            ha="left",
            fontsize=8,
        )

    ax.set_xlabel("Number of PCA dimensions")
    ax.set_ylabel("k-NN overlap")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout(pad=0.6)

    pdf_path = OUT_DIR / "sirna_pca_dependence.pdf"
    fig.savefig(pdf_path)
    plt.close(fig)

    print(f"Saved: {pdf_path}")
    print(f"Saved: {source_data}")


if __name__ == "__main__":
    main()
