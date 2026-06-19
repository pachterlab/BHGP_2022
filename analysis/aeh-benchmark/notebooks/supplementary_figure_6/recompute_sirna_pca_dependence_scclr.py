#!/usr/bin/env python3
"""Recompute Supplementary Note Figure 7 PFlog rows with scclr.

The checked-in Ahlmann-Eltze--Huber downsampling table contains the original
``clr`` rows. Figure 2c replaces those rows with scclr sparse shifted-CLR PCA
using estimated alpha/K; this script does the same for the siRNA PCA-dimension
sweep shown in Supplementary Note Figure 7.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parents[3]
FIGURE_DIR = REPO_DIR.parent / "figures"
OUT = SCRIPT_DIR / "sirna_pca_dependence_clr_alpha_k.tsv"


def load_fig2_helpers():
    helper_path = REPO_DIR / "scripts" / "recompute_fig2abc_scclr.py"
    spec = importlib.util.spec_from_file_location("recompute_fig2abc_scclr", helper_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {helper_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ensure_dirs()
    return module


def main() -> None:
    helpers = load_fig2_helpers()
    dataset = "smartSeq3_siRNA_knockdown"
    pca_dims = [5, 10, 50, 100]
    seeds = range(1, 6)
    knn_k = 50

    matrix = helpers.filter_gene_cell(helpers.load_downsampling_gene_cell(dataset))
    rows = []
    for seed in seeds:
        reduced = helpers.r_downsample_matrix(matrix, dataset, seed)
        gene_keep = np.asarray(reduced.sum(axis=1)).ravel() > 0
        X_full = matrix[gene_keep].T.tocsr()
        X_reduced = reduced[gene_keep].T.tocsr()
        cell_keep = (
            (np.asarray(X_full.sum(axis=1)).ravel() > 0)
            & (np.asarray(X_reduced.sum(axis=1)).ravel() > 0)
        )
        X_full = X_full[cell_keep]
        X_reduced = X_reduced[cell_keep]

        for pca_dim in pca_dims:
            emb_full, alpha_full, K_full = helpers.scclr_pca_cells_genes(X_full, pca_dim, seed)
            emb_reduced, alpha_reduced, K_reduced = helpers.scclr_pca_cells_genes(
                X_reduced, pca_dim, seed
            )
            overlap = helpers.mean_overlap(
                helpers.knn(emb_full, knn_k),
                helpers.knn(emb_reduced, knn_k),
            )
            rows.append(
                {
                    "overlap": overlap,
                    "dataset": dataset,
                    "seed": seed,
                    "pca_dim": pca_dim,
                    "knn": knn_k,
                    "transformation": "clr",
                    "alpha": "TRUE",
                    "alpha_est_full": alpha_full,
                    "alpha_est_reduced": alpha_reduced,
                    "K_full": K_full,
                    "K_reduced": K_reduced,
                }
            )
            print(
                f"{dataset} seed={seed} pca_dim={pca_dim} "
                f"overlap={overlap:.3f} K=({K_full:.3g},{K_reduced:.3g})",
                flush=True,
            )

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)

    # Also write a copy next to the Figure 2 override tables for easier auditing.
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(FIGURE_DIR / OUT.name, sep="\t", index=False)
    print(f"Saved: {OUT}")
    print(f"Saved: {FIGURE_DIR / OUT.name}")


if __name__ == "__main__":
    main()
