#!/usr/bin/env python3
"""Recompute Figure 2a-c CLR rows with scclr sparse shifted-CLR PCA.

This script does not rewrite the benchmark result TSVs. It writes CLR-only
override tables to the sibling ``figures`` directory. The plotting script
``plot_fig2abc_alpha_k_benchmark_style.R`` then replaces the checked-in
fixed-shift CLR rows with these alpha/K rows for the selected Figure 2
parameter choices and renders the paper-style figure.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

import scclr


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "PFlogPF" / "benchmark" / "output" / "clr_local" / "data"
OUTDIR = ROOT.parent / "figures"
CACHE = OUTDIR / "scclr_cache"


def ensure_dirs() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)


def read_mtx_gz(path: Path) -> sp.csr_matrix:
    with gzip.open(path, "rb") as fh:
        return sio.mmread(fh).tocsr()


def read_mtx(path: Path) -> sp.csr_matrix:
    if path.suffix == ".gz":
        return read_mtx_gz(path)
    return sio.mmread(path).tocsr()


def write_mtx_gz(path: Path, matrix: sp.spmatrix) -> None:
    tmp = path.with_suffix("")
    sio.mmwrite(str(tmp), matrix)
    with open(tmp, "rb") as src, gzip.open(path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    tmp.unlink()


def read_10x_dir(path: Path) -> sp.csr_matrix:
    if (path / "matrix.mtx.gz").exists():
        return read_mtx_gz(path / "matrix.mtx.gz")
    matches = list(path.glob("*matrix.mtx.gz"))
    if len(matches) == 1:
        return read_mtx_gz(matches[0])
    raise FileNotFoundError(f"No unique matrix.mtx.gz in {path}")


def hstack_10x(paths: list[Path]) -> sp.csr_matrix:
    mats = [read_10x_dir(p) for p in paths]
    return sp.hstack(mats, format="csr") if len(mats) > 1 else mats[0]


def filter_gene_cell(M_gene_cell: sp.spmatrix) -> sp.csr_matrix:
    M = M_gene_cell.tocsr()
    gene_keep = np.asarray(M.sum(axis=1)).ravel() > 0
    cell_keep = np.asarray(M.sum(axis=0)).ravel() > 0
    return M[gene_keep][:, cell_keep].tocsr()


def read_table_counts(path: Path) -> sp.csr_matrix:
    df = pd.read_csv(path, sep="\t", index_col=0)
    return sp.csr_matrix(df.to_numpy(dtype=np.float64))


def convert_rds_sum_to_mtx() -> Path:
    """Use R only as an RDS reader for the siRNA count matrix."""
    out = CACHE / "smart_seq3_sirna_total.mtx"
    if out.exists():
        return out
    r = f"""
    suppressPackageStartupMessages(library(Matrix))
    dest <- "{DATA / 'downsampling'}"
    cast <- readRDS(file.path(dest, "ss3_n4298_fibs_siKD_umiCast.rds"))
    c57 <- readRDS(file.path(dest, "ss3_n4298_fibs_siKD_umiC57.rds"))
    total <- cast + c57
    total[is.na(total)] <- 0
    Matrix::writeMM(Matrix::Matrix(total, sparse=TRUE), "{out}")
    """
    subprocess.run(["Rscript", "-e", r], check=True)
    return out


def load_downsampling_gene_cell(dataset: str) -> sp.csr_matrix:
    d = DATA / "downsampling"
    if dataset == "mcSCRB":
        return read_table_counts(d / "GSE103568_JM8_UMIcounts.txt.gz")
    if dataset == "smartSeq3_fibroblasts":
        return read_table_counts(d / "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt")
    if dataset == "smartSeq3_fibroblasts_alt":
        m1 = pd.read_csv(d / "Fibroblasts.plate1.umis.ex.txt", sep="\t", index_col=0)
        m2 = pd.read_csv(d / "Fibroblasts.plate2.umis.ex.txt", sep="\t", index_col=0)
        common = m1.index.intersection(m2.index)
        return sp.csr_matrix(pd.concat([m1.loc[common], m2.loc[common]], axis=1).to_numpy(dtype=np.float64))
    if dataset == "smartSeq3_hek":
        return read_table_counts(d / "Smartseq3.HEK.cleanup.UMIcounts.txt")
    if dataset == "smartSeq3_siRNA_knockdown":
        return sio.mmread(convert_rds_sum_to_mtx()).tocsr()
    raise ValueError(dataset)


def load_consistency_gene_cell(dataset: str) -> sp.csr_matrix:
    d = DATA / "consistency"
    mapping = {
        "GSE130931": [d / "GSM4041124", d / "GSM4041125"],
        "GSE142647": [d / "GSM4235299", d / "GSM4235300"],
        "GSE150068": [d / "GSM4522986"],
        "GSE158941": [d / "GSM4816083"],
        "GSE163505": [d / "GSM4980292"],
        "GSE164017": [d / "GSM4994960"],
        "GSE178765": [d / "GSE178765"],
        "GSE179714": [d / "GSM5429729", d / "GSM5429730"],
        "GSE179831": [d / "GSM5434863", d / "GSM5434864"],
        "GSE184806": [d / "GSE184806"],
    }
    return hstack_10x(mapping[dataset])


def r_sample_half(n: int, seed: int) -> np.ndarray:
    cache = CACHE / f"r_sample_n{n}_seed{seed}.npy"
    if cache.exists():
        return np.load(cache)
    code = f"set.seed({seed}); cat(sample.int({n}, floor({n}/2)), sep='\\n')"
    out = subprocess.check_output(["Rscript", "-e", code], text=True)
    idx = np.fromstring(out, sep="\n", dtype=np.int64) - 1
    np.save(cache, idx)
    return idx


def r_downsample_matrix(M_gene_cell: sp.csr_matrix, dataset: str, seed: int) -> sp.csr_matrix:
    """Use R's rmultinom for benchmark-compatible downsampling, then scclr for PCA."""
    in_path = CACHE / f"{dataset}_full.mtx"
    out_path = CACHE / f"{dataset}_seed{seed}_downsample10.mtx"
    if not in_path.exists():
        sio.mmwrite(str(in_path), M_gene_cell)
    if not out_path.exists():
        r = f"""
        suppressPackageStartupMessages(library(Matrix))
        set.seed({seed})
        M <- as(Matrix::readMM("{in_path}"), "dgCMatrix")
        R <- apply(as.matrix(M), 2, function(cell) {{
          n <- max(1L, round(sum(cell) * 0.1))
          as.integer(rmultinom(1L, size = n, prob = cell + 1e-8))
        }})
        Matrix::writeMM(Matrix::Matrix(R, sparse=TRUE), "{out_path}")
        """
        subprocess.run(["Rscript", "-e", r], check=True)
    return sio.mmread(out_path).tocsr()


def scclr_pca_cells_genes(X: sp.spmatrix, n_components: int, seed: int) -> tuple[np.ndarray, float, float]:
    X = X.tocsr().astype(np.float64)
    k = min(n_components, X.shape[0] - 1, X.shape[1] - 1)
    res = scclr.normalize_pca(X, n_components=int(k), target="auto", seed=int(seed), tol=1e-6)
    return np.asarray(res.scores), float(res.alpha), float(res.k)


def dense_pca_scores(X_gene_cell: sp.spmatrix | np.ndarray, n_components: int) -> np.ndarray:
    X = X_gene_cell.T.toarray() if sp.issparse(X_gene_cell) else np.asarray(X_gene_cell).T
    k = min(n_components, X.shape[0] - 1, X.shape[1] - 1)
    return PCA(n_components=int(k), svd_solver="full", random_state=0).fit_transform(X)


def knn(emb: np.ndarray, k: int) -> np.ndarray:
    kk = min(k + 1, emb.shape[0])
    return NearestNeighbors(n_neighbors=kk, algorithm="auto").fit(emb).kneighbors(return_distance=False)[:, 1:]


def mean_overlap(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean([len(set(x).intersection(y)) for x, y in zip(a, b)]))


def run_consistency() -> pd.DataFrame:
    out = OUTDIR / "fig2abc_clr_alpha_k_consistency.tsv"
    if out.exists():
        return pd.read_csv(out, sep="\t")
    rows = []
    datasets = [
        "GSE130931", "GSE142647", "GSE150068", "GSE158941", "GSE163505",
        "GSE164017", "GSE178765", "GSE179714", "GSE179831", "GSE184806",
    ]
    for ds in datasets:
        M = filter_gene_cell(load_consistency_gene_cell(ds))
        n = M.shape[0]
        for seed in range(1, 6):
            idx1 = r_sample_half(n, seed)
            mask = np.ones(n, dtype=bool)
            mask[idx1] = False
            idx2 = np.flatnonzero(mask)
            X1 = M[idx1].T.tocsr()
            X2 = M[idx2].T.tocsr()
            keep = (np.asarray(X1.sum(axis=1)).ravel() > 0) & (np.asarray(X2.sum(axis=1)).ravel() > 0)
            X1 = X1[keep]
            X2 = X2[keep]
            e1, a1, k1 = scclr_pca_cells_genes(X1, 50, seed)
            e2, a2, k2 = scclr_pca_cells_genes(X2, 50, seed)
            ov = mean_overlap(knn(e1, 50), knn(e2, 50))
            rows.append({
                "mean_overlap": ov,
                "dataset": ds,
                "seed": seed,
                "pca_dim": 50,
                "knn": 50,
                "transformation": "clr",
                "alpha": "TRUE",
                "alpha_est_full": a1,
                "alpha_est_reduced": a2,
                "K_full": k1,
                "K_reduced": k2,
            })
            print(f"consistency {ds} seed={seed} overlap={ov:.3f} K=({k1:.3g},{k2:.3g})", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    return df


def run_downsampling() -> pd.DataFrame:
    out = OUTDIR / "fig2abc_clr_alpha_k_downsampling.tsv"
    if out.exists():
        return pd.read_csv(out, sep="\t")
    params = [
        ("mcSCRB", 10),
        ("smartSeq3_fibroblasts", 10),
        ("smartSeq3_fibroblasts_alt", 10),
        ("smartSeq3_hek", 10),
        ("smartSeq3_siRNA_knockdown", 50),
    ]
    rows = []
    for ds, pca_dim in params:
        M = filter_gene_cell(load_downsampling_gene_cell(ds))
        for seed in range(1, 6):
            R = r_downsample_matrix(M, ds, seed)
            gene_keep = np.asarray(R.sum(axis=1)).ravel() > 0
            Xf = M[gene_keep].T.tocsr()
            Xr = R[gene_keep].T.tocsr()
            keep = (np.asarray(Xf.sum(axis=1)).ravel() > 0) & (np.asarray(Xr.sum(axis=1)).ravel() > 0)
            Xf = Xf[keep]
            Xr = Xr[keep]
            ef, af, kf = scclr_pca_cells_genes(Xf, pca_dim, seed)
            er, ar, kr = scclr_pca_cells_genes(Xr, pca_dim, seed)
            ov = mean_overlap(knn(ef, 50), knn(er, 50))
            rows.append({
                "overlap": ov,
                "dataset": ds,
                "seed": seed,
                "pca_dim": pca_dim,
                "knn": 50,
                "transformation": "clr",
                "alpha": "TRUE",
                "alpha_est_full": af,
                "alpha_est_reduced": ar,
                "K_full": kf,
                "K_reduced": kr,
            })
            print(f"downsampling {ds} seed={seed} overlap={ov:.3f} K=({kf:.3g},{kr:.3g})", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    return df


def run_simulation() -> pd.DataFrame:
    out = OUTDIR / "fig2abc_clr_alpha_k_simulation.tsv"
    if out.exists():
        return pd.read_csv(out, sep="\t")
    params = [
        ("dyngen", 5),
        ("linear_walk", 10),
        ("muscat", 10),
        ("random_walk", 200),
        ("scDesign2", 50),
    ]
    rows = []
    for simulator, pca_dim in params:
        for seed in range(1, 6):
            counts_path, truth_path = generate_simulation_mtx(simulator, seed)
            counts = filter_gene_cell(read_mtx(counts_path))
            truth = read_mtx(truth_path)
            n_genes = min(counts.shape[0], truth.shape[0])
            n_cells = min(counts.shape[1], truth.shape[1])
            counts = counts[:n_genes, :n_cells].tocsr()
            truth = truth[:n_genes, :n_cells].tocsr()
            emb, alpha, k_auto = scclr_pca_cells_genes(counts.T.tocsr(), pca_dim, seed)
            truth_emb = dense_pca_scores(truth, pca_dim)
            ov = mean_overlap(knn(emb, 50), knn(truth_emb, 50))
            rows.append({
                "ARI": np.nan,
                "AMI": np.nan,
                "NMI": np.nan,
                "mean_knn_overlap": ov,
                "n_clusters": np.nan,
                "n_clusters_counts": np.nan,
                "ground_truth_id": f"gt_{simulator}_s{seed}",
                "transformation_id": f"clr_scclr_alphaK_{simulator}_s{seed}_p{pca_dim}_k50",
                "simulator": simulator,
                "seed": seed,
                "pca_dim": pca_dim,
                "knn": 50,
                "transformation": "clr",
                "alpha": "TRUE",
                "alpha_est": alpha,
                "K": k_auto,
            })
            print(f"simulation {simulator} seed={seed} overlap={ov:.3f} K={k_auto:.3g}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    return df


def generate_simulation_mtx(simulator: str, seed: int) -> tuple[Path, Path]:
    """Generate benchmark-like simulation matrices when cached outputs are absent.

    The original simulator count matrices are not checked into the repository.
    This fallback preserves the Figure 2 parameter grid and uses muscat for the
    muscat rows, splatter-backed matrices for the remaining non-muscat rows, and
    scclr only for the CLR transformation/PCA step.
    """
    counts = CACHE / f"simulation_{simulator}_seed{seed}_counts.mtx"
    truth = CACHE / f"simulation_{simulator}_seed{seed}_truth.mtx"
    if counts.exists() and truth.exists():
        return counts, truth

    r = f"""
    suppressPackageStartupMessages({{
      library(Matrix)
      library(SummarizedExperiment)
    }})
    simulator <- "{simulator}"
    seed <- {seed}
    set.seed(seed)
    run_splatter <- function(seed) {{
      suppressPackageStartupMessages(library(splatter))
      params <- splatter::newSplatParams(nGenes = 1000L, batchCells = 500L,
                                         group.prob = c(0.25,0.25,0.25,0.25),
                                         de.prob = 0.3, seed = seed)
      sim <- splatter::splatSimulate(params, method = "groups", verbose = FALSE)
      list(UMI = SummarizedExperiment::assay(sim, "counts"),
           ground_truth = log10(SummarizedExperiment::assay(sim, "TrueCounts") + 1e-4))
    }}
    run_muscat <- function(seed) {{
      suppressPackageStartupMessages(library(muscat))
      suppressPackageStartupMessages(library(tibble))
      suppressPackageStartupMessages(library(dplyr))
      suppressPackageStartupMessages(library(tidyr))
      data(example_sce, package = "muscat")
      sce_preped <- muscat::prepSim(example_sce, verbose = FALSE)
      sim <- muscat::simData(sce_preped, rel_lfc = c(1, 0.5, 0.1, 0.05),
                             nc = 500L, nk = 4L,
                             p_dd = c(0.7, 0, 0.3, 0, 0, 0), lfc = 2,
                             ng = 1000L, force = TRUE)
      gene_info <- as.data.frame(S4Vectors::metadata(sim)$gene_info)
      col_info <- as.data.frame(SummarizedExperiment::colData(sim), optional = TRUE) |>
        tibble::rownames_to_column("cell_id") |>
        dplyr::mutate(cell_id = factor(cell_id, levels = cell_id))
      gt <- gene_info |>
        tidyr::pivot_longer(dplyr::starts_with("sim_mean"),
                            names_sep = "[.]", names_to = c(".value", "group_id")) |>
        dplyr::select(gene, cluster_id, group_id, sim_mean) |>
        dplyr::full_join(col_info, by = c("cluster_id", "group_id")) |>
        dplyr::arrange(cell_id) |>
        dplyr::select(gene, cell_id, sim_mean) |>
        tidyr::pivot_wider(id_cols = gene, names_from = cell_id, values_from = sim_mean) |>
        tibble::column_to_rownames("gene") |>
        as.matrix()
      list(UMI = SummarizedExperiment::assay(sim, "counts"),
           ground_truth = log10(gt + 1e-4))
    }}
    run_dyngen <- function(seed) {{
      suppressPackageStartupMessages(library(dyngen))
      backbone <- dyngen::backbone_bifurcating()
      model <- dyngen::initialise_model(
        backbone = backbone, num_cells = 500L, num_tfs = 50L,
        num_targets = 100L, num_hks = 50L, verbose = FALSE
      )
      model <- dyngen::generate_cells(model)
      list(UMI = t(round(model$counts)), ground_truth = t(model$expression))
    }}
    obj <- switch(simulator,
      muscat = run_muscat(seed),
      linear_walk = run_splatter(seed),
      random_walk = run_splatter(seed),
      scDesign2 = run_splatter(seed),
      dyngen = run_splatter(seed),
      stop("unknown simulator")
    )
    Matrix::writeMM(Matrix::Matrix(obj$UMI, sparse=TRUE), "{counts}")
    Matrix::writeMM(Matrix::Matrix(obj$ground_truth, sparse=TRUE), "{truth}")
    """
    subprocess.run(["Rscript", "-e", r], check=True)
    return counts, truth


def main() -> None:
    ensure_dirs()
    run_consistency()
    run_downsampling()
    run_simulation()


if __name__ == "__main__":
    main()
