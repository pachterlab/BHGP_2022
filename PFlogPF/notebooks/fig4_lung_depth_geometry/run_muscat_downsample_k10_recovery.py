#!/usr/bin/env python3
"""Compute the muscat 12x downsampling k=10 neighborhood-recovery result."""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
DEFAULT_COUNTS = HERE / "data" / "muscat_seed1_counts.mtx"
DEFAULT_LABELS = HERE / "data" / "muscat_seed1_truth_cluster_labels.tsv"
DEFAULT_OUTPUT = HERE / "tables" / "muscat_downsample12_k10_recovery.tsv"


def import_scclr(extra_path: Path | None = None):
    if extra_path is not None:
        sys.path.insert(0, str(extra_path))
    import scclr

    return scclr


def import_scvi():
    try:
        import scvi
        import torch
    except ImportError as exc:
        raise SystemExit(
            "This script requires scvi-tools and torch for the scVI/scANVI rows. "
            "Install scvi-tools or run in the environment used for manuscript analysis."
        ) from exc
    return scvi, torch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--scclr-path", type=Path, default=None, help="Optional path to a local scclr/python directory.")
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--downsample-probability", type=float, default=1 / 12)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--scvi-max-epochs", type=int, default=100)
    parser.add_argument("--scanvi-max-epochs", type=int, default=60)
    return parser.parse_args()


def set_seed(scvi, torch, seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    scvi.settings.seed = seed


def downsample_csr(X: sp.csr_matrix, p: float, seed: int) -> sp.csr_matrix:
    rng = np.random.default_rng(seed)
    Y = X.copy().astype(np.int64)
    Y.data = rng.binomial(Y.data, p)
    Y.eliminate_zeros()
    return Y.tocsr()


def log_cp10k(X: sp.csr_matrix) -> np.ndarray:
    depth = np.asarray(X.sum(axis=1)).ravel()
    scale = np.divide(10000.0, depth, out=np.zeros_like(depth, dtype=float), where=depth > 0)
    Y = X.multiply(scale[:, None]).tocsr()
    Y.data = np.log1p(Y.data)
    return Y.toarray()


def pca_scores(Y: np.ndarray, n_components: int = 30) -> np.ndarray:
    n_components = min(n_components, Y.shape[0] - 1, Y.shape[1])
    return PCA(n_components=n_components, svd_solver="full").fit_transform(Y - Y.mean(axis=0, keepdims=True))


def full_knn_recovery_from_downsampled(emb: np.ndarray, k: int) -> float:
    n = emb.shape[0] // 2
    high = StandardScaler().fit_transform(emb[:n])
    low = StandardScaler().fit_transform(emb[n:])
    reference = NearestNeighbors(n_neighbors=k + 1).fit(high).kneighbors(high, return_distance=False)[:, 1:]
    predicted_with_twin = NearestNeighbors(n_neighbors=k + 1).fit(high).kneighbors(low, return_distance=False)

    vals = []
    for i in range(n):
        predicted = [j for j in predicted_with_twin[i] if j != i][:k]
        vals.append(len(set(reference[i]) & set(predicted)) / k)
    return float(np.mean(vals))


def fit_scvi_scanvi(scvi, joint: sp.csr_matrix, batch: np.ndarray, labels: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    adata = ad.AnnData(joint.astype(np.float32))
    adata.obs["batch"] = batch
    adata.obs["truth_label"] = labels
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="truth_label")
    vae = scvi.model.SCVI(adata, n_latent=30, n_layers=2, n_hidden=128, gene_likelihood="nb")
    vae.train(max_epochs=args.scvi_max_epochs, check_val_every_n_epoch=20, enable_progress_bar=False)
    scvi_latent = vae.get_latent_representation()

    scanvi = scvi.model.SCANVI.from_scvi_model(vae, labels_key="truth_label", unlabeled_category="Unknown")
    scanvi.train(max_epochs=args.scanvi_max_epochs, check_val_every_n_epoch=20, enable_progress_bar=False)
    scanvi_latent = scanvi.get_latent_representation()
    return scvi_latent, scanvi_latent


def main() -> None:
    args = parse_args()
    scclr = import_scclr(args.scclr_path)
    scvi, torch = import_scvi()
    set_seed(scvi, torch, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    counts = mmread(args.counts).T.tocsr().astype(np.int64)
    low = downsample_csr(counts, p=args.downsample_probability, seed=args.seed)
    joint = sp.vstack([counts, low], format="csr")
    n = counts.shape[0]

    # scANVI requires labels. The default labels are coarse k-means clusters
    # computed from the muscat truth matrix used in the simulator benchmark.
    if args.labels.exists():
        labels = pd.read_csv(args.labels, sep="\t")["truth_cluster"].to_numpy()
    else:
        full_pca = pca_scores(log_cp10k(counts), n_components=20)
        labels = KMeans(n_clusters=4, n_init=50, random_state=args.seed).fit_predict(StandardScaler().fit_transform(full_pca))
    joint_labels = np.array([f"cluster_{x}" for x in np.r_[labels, labels]])
    batch = np.array(["full"] * n + ["downsample12x"] * n)

    embeddings: dict[str, np.ndarray] = {
        "logCP10K PCA": pca_scores(log_cp10k(joint), n_components=30),
        "PFlog PCA": scclr.normalize_pca(joint, n_components=30, target="auto", seed=args.seed).scores,
    }
    embeddings["scVI latent"], embeddings["scANVI latent"] = fit_scvi_scanvi(scvi, joint, batch, joint_labels, args)

    rows = []
    for method, emb in embeddings.items():
        rows.append(
            {
                "method": method,
                "k": args.k,
                "full_knn_recovery_from_downsampled": full_knn_recovery_from_downsampled(emb, k=args.k),
                "n_cells_original": n,
                "n_genes": counts.shape[1],
                "full_mean_depth": float(np.asarray(counts.sum(axis=1)).mean()),
                "downsampled_mean_depth": float(np.asarray(low.sum(axis=1)).mean()),
                "downsample_probability": args.downsample_probability,
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(args.output, sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"wrote: {args.output}")


if __name__ == "__main__":
    main()
