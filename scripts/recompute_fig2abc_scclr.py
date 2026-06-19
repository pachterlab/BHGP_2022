#!/usr/bin/env python3
"""Recompute adapted AEH benchmark PFlog rows with scclr sparse shifted-CLR PCA.

This script does not rewrite the benchmark result TSVs. It writes CLR-only
override tables to the sibling ``figures`` directory. The plotting script
``plot_fig2abc_alpha_k_benchmark_style.R`` then replaces the checked-in
fixed-shift CLR rows with these alpha/K rows for the selected Supplementary Note
parameter choices and renders the paper-style figure.
"""

from __future__ import annotations

import gzip
import os
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
DATA = ROOT / "analysis" / "aeh-benchmark" / "benchmark" / "output" / "clr_local" / "data"
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
        ("muscat", 10),
        ("scDesign2", 50),
    ]
    rows = []
    for simulator, pca_dim in params:
        for seed in range(1, 6):
            counts_path, truth_knn_path = generate_simulation_mtx(simulator, seed)
            counts_raw = read_mtx(counts_path).tocsr()
            gene_keep = np.asarray(counts_raw.sum(axis=1)).ravel() > 0
            cell_keep = np.asarray(counts_raw.sum(axis=0)).ravel() > 0
            counts = counts_raw[gene_keep][:, cell_keep].tocsr()
            emb, alpha, k_auto = scclr_pca_cells_genes(counts.T.tocsr(), pca_dim, seed)
            # The original benchmark constructs the ground-truth kNN directly in
            # simulator truth space with BiocNeighbors::findAnnoy(..., k=50).
            truth_knn = np.loadtxt(truth_knn_path, delimiter="\t", dtype=np.int64) - 1
            truth_knn = np.atleast_2d(truth_knn)[cell_keep]
            ov = mean_overlap(knn(emb, 50), truth_knn)
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
    """Generate simulation matrices with all original benchmark simulator families.

    The Baron-reference linear/random-walk simulators are capped locally at the
    same 1000-gene/5000-cell scale used by muscat and dyngen; otherwise their
    dense truth space makes local regeneration impractically slow.
    """
    counts = CACHE / f"simulation_allv2_{simulator}_seed{seed}_counts.mtx"
    truth_knn = CACHE / f"simulation_allv2_{simulator}_seed{seed}_truth_knn.tsv"
    if counts.exists() and truth_knn.exists():
        return counts, truth_knn

    legacy_truth = CACHE / f"simulation_allv2_{simulator}_seed{seed}_truth.mtx"
    if counts.exists() and legacy_truth.exists() and not truth_knn.exists():
        r_convert = f"""
        suppressPackageStartupMessages({{
          library(BiocNeighbors)
          library(Matrix)
        }})
        truth <- Matrix::readMM("{legacy_truth}")
        truth_knn <- BiocNeighbors::findKNN(
          as.matrix(t(truth)),
          k = 50L,
          BNPARAM = BiocNeighbors::AnnoyParam(),
          get.distance = FALSE
        )$index
        write.table(truth_knn, file = "{truth_knn}", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
        r_convert_script = CACHE / f"convert_exact_{simulator}_seed{seed}_truth_knn.R"
        r_convert_script.write_text(r_convert)
        subprocess.run(["Rscript", "--vanilla", str(r_convert_script)], check=True)
        return counts, truth_knn

    benchmark_data = DATA
    gse130931_dirs = [
        benchmark_data / "consistency" / "GSM4041124",
        benchmark_data / "consistency" / "GSM4041125",
    ]
    r = f"""
    suppressPackageStartupMessages({{
      library(BiocNeighbors)
      library(Matrix)
      library(MatrixGenerics)
      library(SingleCellExperiment)
      library(SummarizedExperiment)
    }})
    simulator <- "{simulator}"
    seed <- {seed}
    set.seed(seed)

    require_pkg <- function(pkg) {{
      if (!requireNamespace(pkg, quietly = TRUE)) {{
        stop("Required R package not installed for exact simulator rerun: ", pkg, call. = FALSE)
      }}
    }}

    run_muscat <- function(seed) {{
      require_pkg("muscat")
      require_pkg("tibble")
      require_pkg("dplyr")
      require_pkg("tidyr")
      suppressPackageStartupMessages(library(muscat))
      suppressPackageStartupMessages(library(tibble))
      suppressPackageStartupMessages(library(dplyr))
      suppressPackageStartupMessages(library(tidyr))
      data(example_sce, package = "muscat")
      sce_preped <- muscat::prepSim(example_sce, verbose = FALSE)
      sim <- muscat::simData(sce_preped, rel_lfc = c(1, 0.5, 0.1, 0.05),
                             nc = 5000L, nk = 4L,
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
      require_pkg("dyngen")
      suppressPackageStartupMessages(library(dyngen))
      options(Ncpus = 1L)
      options(dyngen_download_cache_dir = "{CACHE / 'dyngen_download_cache'}")
      dir.create(getOption("dyngen_download_cache_dir"), recursive = TRUE, showWarnings = FALSE)
      num_cells <- 5000L
      num_features <- 1000L
      backbone <- dyngen::backbone_consecutive_bifurcating()
      num_tfs <- nrow(backbone$module_info)
      num_targets <- round((num_features - num_tfs) / 2)
      num_hks <- num_features - num_targets - num_tfs
      config <- dyngen::initialise_model(
        backbone = backbone,
        num_tfs = num_tfs,
        num_targets = num_targets,
        num_hks = num_hks,
        num_cells = num_cells
      )
      res <- dyngen::generate_dataset(config, format = "sce")
      sim <- res$dataset
      UMI <- SummarizedExperiment::assay(sim, "counts")
      expressed_cells <- MatrixGenerics::colSums2(UMI) > 10
      expressed_genes <- MatrixGenerics::rowSums2(UMI) > 0
      sim <- sim[expressed_genes, expressed_cells]
      list(UMI = SummarizedExperiment::assay(sim, "counts"),
           ground_truth = t(SingleCellExperiment::reducedDim(sim, "MDS")))
    }}

    load_baron_counts <- function() {{
      require_pkg("scRNAseq")
      require_pkg("scuttle")
      require_pkg("scran")
      sce <- scRNAseq::BaronPancreasData("human")
      sce <- scuttle::logNormCounts(sce)
      SummarizedExperiment::colData(sce)$cluster_id <- scran::quickCluster(sce)
      reference_data_counts <- SummarizedExperiment::assay(sce)
      reference_data_counts <- reference_data_counts[, MatrixGenerics::colSums2(reference_data_counts) > 0]
      reference_data_counts <- reference_data_counts[MatrixGenerics::rowSums2(reference_data_counts) > 0, ]
      reference_data_counts <- reference_data_counts[
        seq_len(min(1000L, nrow(reference_data_counts))),
        seq_len(min(5000L, ncol(reference_data_counts)))
      ]
      reference_data_counts
    }}

    sanity_counts_from_delta <- function(reference_data_counts, delta_true) {{
      n_genes <- nrow(reference_data_counts)
      n_cells <- ncol(reference_data_counts)
      N_c <- MatrixGenerics::colSums2(reference_data_counts)
      mu_tilde_g <- log(MatrixGenerics::rowSums2(reference_data_counts) / sum(reference_data_counts))
      sig2_g <- rexp(n_genes, rate = 1 / 2)
      lambda <- sqrt(sig2_g / matrixStats::rowVars(delta_true))
      delta_true <- (delta_true - MatrixGenerics::rowMeans2(delta_true)) * lambda
      mu_g <- mu_tilde_g - sig2_g / 2
      UMI <- matrix(
        rnbinom(n_genes * n_cells, mu = t(t(exp(mu_g + delta_true)) * N_c), size = 1 / 0.01),
        n_genes,
        n_cells
      )
      list(UMI = UMI, ground_truth = delta_true)
    }}

    run_linear_walk <- function(seed) {{
      reference_data_counts <- load_baron_counts()
      n_genes <- nrow(reference_data_counts)
      n_cells <- ncol(reference_data_counts)
      delta_true <- matrix(NA_real_, n_genes, n_cells)
      parents <- rep(NA_integer_, n_cells)
      branch_length <- 1200L
      branch_idx <- 0L
      start_point <- NULL
      end_point <- NULL
      for (idx in seq_len(n_cells) - 1L) {{
        if (idx == 0L) {{
          parents[idx + 1L] <- -1L
          start_point <- rnorm(n_genes, mean = 0, sd = 1)
          end_point <- rnorm(n_genes, mean = 0, sd = 1)
        }} else if (idx %% branch_length == 0L) {{
          start_id <- sample.int(idx, size = 1L)
          parents[idx + 1L] <- start_id - 1L
          branch_idx <- 0L
          start_point <- rnorm(n_genes, mean = delta_true[, parents[idx + 1L]], sd = 0.1)
          end_point <- rnorm(n_genes, mean = 0, sd = 5)
        }} else {{
          branch_idx <- branch_idx + 1L
          parents[idx + 1L] <- idx - 1L
        }}
        delta_true[, idx + 1L] <- rnorm(
          n_genes,
          mean = start_point + (end_point - start_point) * branch_idx / branch_length,
          sd = 0.1
        )
      }}
      sanity_counts_from_delta(reference_data_counts, delta_true)
    }}

    run_random_walk <- function(seed) {{
      reference_data_counts <- load_baron_counts()
      n_genes <- nrow(reference_data_counts)
      n_cells <- ncol(reference_data_counts)
      delta_true <- matrix(NA_real_, n_genes, n_cells)
      parents <- rep(NA_integer_, n_cells)
      branch_length <- 13L
      for (idx in seq_len(n_cells)) {{
        if (idx == 1L) {{
          delta <- rnorm(n_genes, mean = 0, sd = 1)
          parents[idx] <- 0L
        }} else if (idx %% branch_length == 0L) {{
          parents[idx] <- sample.int(idx - 1L, size = 1L)
          delta <- rnorm(n_genes, mean = delta_true[, parents[idx]], sd = 1)
        }} else {{
          parents[idx] <- idx - 1L
          delta <- rnorm(n_genes, mean = delta_true[, parents[idx]], sd = 1)
        }}
        delta_true[, idx] <- delta
      }}
      sanity_counts_from_delta(reference_data_counts, delta_true)
    }}

    run_scdesign2 <- function(seed) {{
      require_pkg("DropletUtils")
      require_pkg("scran")
      require_pkg("scDesign2")
      dirs <- c("{gse130931_dirs[0]}", "{gse130931_dirs[1]}")
      sce <- DropletUtils::read10xCounts(dirs)
      SummarizedExperiment::colData(sce)$cluster_id <- scran::quickCluster(sce, min.size = 20)
      sce <- sce[, MatrixGenerics::colSums2(SummarizedExperiment::assay(sce)) > 0]
      sce <- sce[MatrixGenerics::rowSums2(SummarizedExperiment::assay(sce)) > 0, ]
      sce <- scuttle::logNormCounts(sce)
      mat <- SummarizedExperiment::assay(sce, "counts")
      colnames(mat) <- as.character(SummarizedExperiment::colData(sce)$cluster_id)
      cluster_names <- unique(as.character(SummarizedExperiment::colData(sce)$cluster_id))
      fit <- scDesign2::fit_model_scDesign2(
        mat,
        cell_type_sel = cluster_names,
        sim_method = "copula",
        marginal = "nb"
      )
      UMI <- scDesign2::simulate_count_scDesign2(
        fit,
        n_cell_new = ncol(mat),
        cell_type_prop = table(colnames(mat))[cluster_names] / ncol(mat),
        sim_method = "copula"
      )
      ground_truth <- do.call(cbind, lapply(fit[cluster_names], function(fi) {{
        gt <- matrix(NA_real_, nrow = nrow(mat), ncol = fi$n_cell)
        gt[fi$gene_sel1, ] <- fi$marginal_param1[, 3]
        gt[fi$gene_sel2, ] <- fi$marginal_param2[, 3]
        gt[fi$gene_sel3, ] <- 1e-8
        gt
      }}))
      list(UMI = UMI, ground_truth = log10(ground_truth))
    }}

    obj <- switch(simulator,
      muscat = run_muscat(seed),
      dyngen = run_dyngen(seed),
      linear_walk = run_linear_walk(seed),
      random_walk = run_random_walk(seed),
      scDesign2 = run_scdesign2(seed),
      stop("unknown simulator")
    )
    Matrix::writeMM(Matrix::Matrix(obj$UMI, sparse=TRUE), "{counts}")
    truth_knn <- BiocNeighbors::findKNN(
      t(obj$ground_truth),
      k = 50L,
      BNPARAM = BiocNeighbors::AnnoyParam(),
      get.distance = FALSE
    )$index
    write.table(truth_knn, file = "{truth_knn}", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
    env = os.environ.copy()
    env.setdefault("XDG_CACHE_HOME", str(CACHE / "xdg_cache"))
    env.setdefault("R_USER_CACHE_DIR", str(CACHE / "r_cache"))
    r_script = CACHE / f"generate_exact_{simulator}_seed{seed}.R"
    r_script.write_text(r)
    subprocess.run(["Rscript", "--vanilla", str(r_script)], check=True, env=env)
    return counts, truth_knn


def main() -> None:
    ensure_dirs()
    run_consistency()
    run_downsampling()
    run_simulation()


if __name__ == "__main__":
    main()
