#!/usr/bin/env python3
"""Compute the Fig. 1a depth negative-control balance metric.

The analysis compares high-depth cells with the immediately lower-depth cells
within the same cell type and age group. Within each biological replicate, the
same number of cells is selected for both depth strata, so the comparison has
identical replicate composition. Genes are abundance matched by binning on raw
mean count and repeatedly forming random gene pairs within bins. For each
normalization, log-ratio balances (z_g - z_h) / sqrt(2) are tested for a
depth-stratum effect after subtracting replicate means.

Example:
    python scripts/metrics_celltype_balance.py \
        --metadata mouse_lung_single_cell_metadata.csv.gz \
        --rdata mouse_lung_single_cell.RData \
        --out-dir analysis/figures
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread, mmwrite

import scclr


METHODS = [
    "raw",
    "PF",
    "sqrt",
    "log1p",
    "log1pCP10k",
    "log1pCPM",
    "sctransform",
    "log1pPF",
    "PFlog",
]


def row_depth(X):
    return np.asarray(X.sum(axis=1)).ravel()


def scale_rows(X, target: float):
    depth = row_depth(X)
    return X.multiply(target / depth[:, None]).tocsr()


def log1p_scaled(X, target: float) -> np.ndarray:
    return scale_rows(X, target).log1p().toarray()


def load_counts(args, out_dir: Path):
    if args.counts_mtx is not None:
        return mmread(args.counts_mtx).tocsr().astype(float)
    if args.rdata is None:
        raise SystemExit("provide either --counts-mtx or --rdata")

    mtx = out_dir / "angelidis_raw_counts_cells_genes.mtx"
    if not mtx.exists():
        code = (
            f'load("{args.rdata}"); '
            f'Matrix::writeMM(Matrix::t(raw_counts), "{mtx}")'
        )
        subprocess.run(["Rscript", "--vanilla", "-e", code], check=True)
    return mmread(mtx).tocsr().astype(float)


def sctransform_residuals(X, out_dir: Path, prefix: str) -> np.ndarray:
    cache = out_dir / f"{prefix}_sctransform_v2.bin"
    shape_file = out_dir / f"{prefix}_sctransform_v2.shape"
    if cache.exists() and shape_file.exists():
        n_cells, n_genes = map(int, shape_file.read_text().strip().split(","))
        return np.fromfile(cache, dtype=np.float32).reshape(n_cells, n_genes)

    mtx = out_dir / f"{prefix}_selected_genes_cells.mtx"
    mmwrite(mtx, X.T.tocsr())
    r_script = out_dir / f"{prefix}_run_sctransform_v2.R"
    r_script.write_text(
        f"""
suppressMessages({{ library(Matrix); library(sctransform) }})
M <- as(readMM("{mtx}"), "CsparseMatrix")
rownames(M) <- paste0("g", seq_len(nrow(M)))
colnames(M) <- paste0("c", seq_len(ncol(M)))
keep_g <- Matrix::rowSums(M) > 0
keep_c <- Matrix::colSums(M) > 0
v <- vst(M[keep_g, keep_c, drop = FALSE], vst.flavor = "v2", verbosity = 0)
y <- v$y
full <- matrix(0.0, nrow = nrow(M), ncol = ncol(M),
               dimnames = list(rownames(M), colnames(M)))
full[rownames(y), colnames(y)] <- y
con <- file("{cache}", "wb")
writeBin(as.vector(full), con, size = 4)
close(con)
writeLines(sprintf("%d,%d", ncol(M), nrow(M)), "{shape_file}")
""".strip()
        + "\n"
    )
    subprocess.run(["Rscript", "--vanilla", str(r_script)], check=True)
    n_cells, n_genes = map(int, shape_file.read_text().strip().split(","))
    return np.fromfile(cache, dtype=np.float32).reshape(n_cells, n_genes)


def residualize_by_group(Y: np.ndarray, group: np.ndarray) -> np.ndarray:
    out = Y.astype(float, copy=True)
    for value in np.unique(group):
        idx = group == value
        out[idx] -= out[idx].mean(axis=0, keepdims=True)
    return out


def depth_test(features: np.ndarray, depth_label: np.ndarray, group: np.ndarray) -> np.ndarray:
    x = residualize_by_group(depth_label.astype(float)[:, None], group).ravel()
    Y = residualize_by_group(features, group)
    xx = float(x @ x)
    df = features.shape[0] - np.unique(group).size - 1
    beta = (x[:, None] * Y).sum(axis=0) / xx
    resid = Y - x[:, None] * beta[None, :]
    sse = (resid * resid).sum(axis=0)
    se = np.sqrt((sse / df) / xx)
    t = np.divide(beta, se, out=np.zeros_like(beta), where=se > 0)
    return 2.0 * stats.t.sf(np.abs(t), df)


def make_pairings(gene_idx: np.ndarray, mean_count: np.ndarray, bins: int, reps: int, seed: int):
    rng = np.random.default_rng(seed)
    sorted_genes = gene_idx[np.argsort(mean_count[gene_idx])]
    abundance_bins = np.array_split(sorted_genes, bins)
    pairings = []
    for _ in range(reps):
        pairs = []
        for b in abundance_bins:
            b = b.copy()
            rng.shuffle(b)
            if b.size % 2:
                b = b[:-1]
            if b.size:
                pairs.append(b.reshape(-1, 2))
        pairings.append(np.vstack(pairs))
    return pairings


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True, help="cell metadata table")
    ap.add_argument("--counts-mtx", default=None, help="cells x genes raw-count Matrix Market file")
    ap.add_argument("--rdata", default=None, help="RData file containing raw_counts genes x cells")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="angelidis_type2_3m_many_balance_de")
    ap.add_argument("--celltype", default="Type_2_pneumocytes")
    ap.add_argument("--celltype-col", default="celltype")
    ap.add_argument("--group-col", default="grouping")
    ap.add_argument("--group-value", default="3m")
    ap.add_argument("--replicate-col", default="identifier")
    ap.add_argument("--n-per-replicate", type=int, default=50)
    ap.add_argument("--detection-threshold", type=float, default=0.25)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--n-reps", type=int, default=200)
    ap.add_argument("--bins", type=int, default=12)
    ap.add_argument("--seed", type=int, default=7)
    return ap.parse_args()


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.metadata)
    counts = load_counts(args, out_dir)
    selected_celltype = (
        meta[args.celltype_col].eq(args.celltype)
        & meta[args.group_col].eq(args.group_value)
    ).to_numpy()
    meta_sub = meta.loc[selected_celltype].copy().reset_index(drop=True)
    X_sub = counts[selected_celltype].tocsr()
    meta_sub["raw_depth"] = row_depth(X_sub)

    low_idx: list[int] = []
    high_idx: list[int] = []
    for _, sub in meta_sub.groupby(args.replicate_col, sort=True):
        order = sub.sort_values("raw_depth", ascending=True).index.to_numpy()
        if order.size < 2 * args.n_per_replicate:
            continue
        low_idx.extend(order[-2 * args.n_per_replicate : -args.n_per_replicate].tolist())
        high_idx.extend(order[-args.n_per_replicate :].tolist())

    low_idx = np.asarray(low_idx)
    high_idx = np.asarray(high_idx)
    selected = np.r_[low_idx, high_idx]
    X = X_sub[selected].tocsr()
    selected_meta = meta_sub.iloc[selected].copy().reset_index(drop=True)
    replicate = selected_meta[args.replicate_col].astype(str).to_numpy()
    depth_label = np.r_[np.zeros(low_idx.size), np.ones(high_idx.size)]
    selected_depth = row_depth(X)
    mean_depth = float(selected_depth.mean())

    methods = {
        "raw": X.toarray(),
        "PF": scale_rows(X, mean_depth).toarray(),
        "sqrt": np.sqrt(X.toarray()),
        "log1p": X.log1p().toarray(),
        "log1pCP10k": log1p_scaled(X, 10_000.0),
        "log1pCPM": log1p_scaled(X, 1_000_000.0),
        "log1pPF": log1p_scaled(X, mean_depth),
    }
    methods["sctransform"] = sctransform_residuals(X, out_dir, args.prefix)
    pflog = scclr.normalize(X, target="auto", center=True)
    methods["PFlog"] = pflog.to_dense()

    n_low = low_idx.size
    detected = (
        (np.asarray((X[:n_low] > 0).sum(axis=0)).ravel() > args.detection_threshold * n_low)
        & (np.asarray((X[n_low:] > 0).sum(axis=0)).ravel() > args.detection_threshold * n_low)
    )
    mean_count = np.asarray(X.mean(axis=0)).ravel()
    gene_idx = np.where(detected)[0]
    pairings = make_pairings(gene_idx, mean_count, args.bins, args.n_reps, args.seed)

    rows = []
    for rep, pairs in enumerate(pairings, start=1):
        for method in METHODS:
            values = methods[method]
            balances = (values[:, pairs[:, 0]] - values[:, pairs[:, 1]]) / np.sqrt(2.0)
            p = depth_test(balances, depth_label, replicate)
            padj = np.minimum(1.0, p * pairs.shape[0])
            rows.append(
                {
                    "replicate": rep,
                    "method": method,
                    "n_balances": int(pairs.shape[0]),
                    "significant_balances": int((padj < args.alpha).sum()),
                    "fraction_significant": float((padj < args.alpha).mean()),
                    "median_neglog10_padj": float(np.median(-np.log10(np.maximum(padj, 1e-300)))),
                    "p90_neglog10_padj": float(np.percentile(-np.log10(np.maximum(padj, 1e-300)), 90)),
                }
            )

    reps = pd.DataFrame(rows)
    reps.to_csv(out_dir / f"{args.prefix}_replicates.tsv", sep="\t", index=False)
    summary = (
        reps.groupby("method", as_index=False)
        .agg(
            mean_significant_balances=("significant_balances", "mean"),
            sd_significant_balances=("significant_balances", "std"),
            median_significant_balances=("significant_balances", "median"),
            mean_fraction_significant=("fraction_significant", "mean"),
            median_p90_neglog10_padj=("p90_neglog10_padj", "median"),
        )
        .sort_values("mean_significant_balances")
    )
    summary.to_csv(out_dir / f"{args.prefix}_summary.tsv", sep="\t", index=False)

    info = {
        "comparison": f"{args.celltype}, {args.group_value}, replicate-balanced high-depth vs next-depth cells",
        "endpoint": "repeated abundance-matched log-ratio balance DE",
        "n_replicates": args.n_reps,
        "bins": args.bins,
        "detection_threshold": args.detection_threshold,
        "n_genes": int(gene_idx.size),
        "n_balances_per_replicate": int(pairings[0].shape[0]),
        "alpha": args.alpha,
        "depth": {
            "low_mean": float(np.mean(selected_depth[:n_low])),
            "high_mean": float(np.mean(selected_depth[n_low:])),
            "low_median": float(np.median(selected_depth[:n_low])),
            "high_median": float(np.median(selected_depth[n_low:])),
        },
        "normalization": {
            "mean_depth": mean_depth,
            "scclr_alpha": float(pflog.alpha),
            "scclr_k": float(pflog.k),
        },
    }
    (out_dir / f"{args.prefix}.json").write_text(json.dumps(info, indent=2) + "\n")


if __name__ == "__main__":
    main()
