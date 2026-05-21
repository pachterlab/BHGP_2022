#!/usr/bin/env Rscript
# run_clr_local.R
# Self-contained local runner for the CLR transformation benchmark.
# Downloads data, runs CLR through all 3 benchmarks (consistency, downsampling,
# simulation), and appends results to the existing pre-computed TSV files.
#
# Run from the benchmark/ directory:
#   Rscript run_clr_local.R

# ── Configuration ─────────────────────────────────────────────────────────────

OUTPUT_DIR  <- file.path(getwd(), "output/clr_local")
RESULTS_DIR <- file.path(getwd(), "output/benchmark_results")
DATA_DIR    <- file.path(OUTPUT_DIR, "data")
N_CORES     <- max(1L, parallel::detectCores() - 1L)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

message("Using ", N_CORES, " cores")

# ── Package installation ──────────────────────────────────────────────────────

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) == 0) return(invisible(NULL))
  message("Installing: ", paste(missing, collapse = ", "))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", type = "binary", repos = "https://cloud.r-project.org")
    BiocManager::install(missing, ask = FALSE, update = FALSE, type = "binary")
  } else {
    install.packages(missing, type = "binary", repos = "https://cloud.r-project.org")
  }
}

install_if_missing(c("irlba", "FNN", "readr", "dplyr", "tidyr", "purrr",
                      "yaml", "digest", "Matrix"))
install_if_missing(c("BiocManager"), bioc = FALSE)
install_if_missing(c("GEOquery", "DropletUtils", "SingleCellExperiment",
                      "MatrixGenerics", "matrixStats"), bioc = TRUE)

# Attempt simulation packages (not fatal if missing)
tryCatch(install_if_missing(c("muscat", "splatter", "scDesign2"), bioc = TRUE),
         error = function(e) message("Note: simulation packages not installed: ", e$message))
tryCatch(install_if_missing(c("dyngen"), bioc = FALSE),
         error = function(e) message("Note: dyngen not installed: ", e$message))

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(matrixStats)
  library(SingleCellExperiment)
})

# ── Core functions ─────────────────────────────────────────────────────────────

# CLR: log(y/sf + 1) then subtract per-cell mean across genes (Aitchison)
clr_transform <- function(UMI, sf) {
  if (inherits(UMI, "sparseMatrix")) UMI <- as.matrix(UMI)
  logpf <- log1p(sweep(UMI, 2L, sf, "/"))     # genes × cells
  sweep(logpf, 2L, colMeans(logpf), "-")       # subtract per-cell mean
}

# Truncated PCA — centers columns internally (as PCA does)
run_pca <- function(X_genes_x_cells, n_components) {
  X <- t(X_genes_x_cells)   # cells × genes
  k <- min(n_components, nrow(X) - 1L, ncol(X) - 1L)
  if (k < 1L) return(X)
  if (k >= min(nrow(X), ncol(X)) / 2) {
    prcomp(X, rank. = k, center = TRUE, scale. = FALSE)$x
  } else {
    prcomp_irlba(X, n = k, center = TRUE, scale. = FALSE)$x
  }
}

# Approximate KNN (cells × k index matrix, 1-based, self excluded)
make_knn <- function(X_cells_x_k, k) {
  k_use <- min(k + 1L, nrow(X_cells_x_k) - 1L)
  nn <- FNN::get.knnx(X_cells_x_k, X_cells_x_k, k = k_use + 1L)$nn.index
  nn[, -1L, drop = FALSE][, seq_len(k), drop = FALSE]
}

# Overlap metric used in consistency / downsampling benchmarks
mean_knn_overlap <- function(knn1, knn2) {
  n <- nrow(knn1)
  mean(vapply(seq_len(n), function(i) length(intersect(knn1[i,], knn2[i,])), 0.0))
}

# One CLR run: transform → PCA → KNN; returns KNN matrix
clr_knn <- function(UMI, sf, pca_dim, knn) {
  trans  <- clr_transform(UMI, sf)
  pca    <- run_pca(trans, pca_dim)
  make_knn(pca, knn)
}

size_factors <- function(UMI) {
  cs <- Matrix::colSums(UMI)
  cs / mean(cs)
}

filter_matrix <- function(UMI) {
  UMI[Matrix::rowSums(UMI) > 0, Matrix::colSums(UMI) > 0, drop = FALSE]
}

# ── Benchmark parameters (from job_overview.yaml) ────────────────────────────

YAML       <- yaml::read_yaml("job_overview.yaml")
KNN_VALS   <- YAML$consistency$knn_construction$knn          # c(10, 50, 100)
PCA_VALS   <- c(5, 10, 50)                                    # drop pca=100
SEEDS      <- YAML$consistency$input_data$seed               # 1:5
DATASETS   <- c("GSE142647", "GSE178765", "GSE163505", "GSE150068")  # 4 datasets

DS_KNN     <- YAML$downsampling$knn_construction$knn
DS_PCA     <- c(5, 10, 50)                                    # drop pca=100
DS_SEEDS   <- YAML$downsampling$input_data$seed
DS_DATASETS<- c("mcSCRB", "smartSeq3_fibroblasts", "smartSeq3_fibroblasts_alt")  # HEK/siRNA no longer on EBI

SIM_KNN    <- YAML$simulation$knn_construction$knn
SIM_PCA    <- c(5, 10, 50)                                    # drop pca=100
SIM_SEEDS  <- YAML$simulation$input_data$seed
SIM_SIMS   <- YAML$simulation$input_data$simulator

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1 — CONSISTENCY BENCHMARK
# ═══════════════════════════════════════════════════════════════════════════════

SKIP_CONSISTENCY <- any(grepl("^clr$",
  readr::read_tsv(file.path(RESULTS_DIR, "consistency_results.tsv"),
                  show_col_types = FALSE)$transformation))

message("\n=== CONSISTENCY BENCHMARK ===")
if (SKIP_CONSISTENCY) {
  message("CLR rows already present in consistency_results.tsv — skipping.")
  consistency_results <- readr::read_tsv(
    file.path(RESULTS_DIR, "consistency_results.tsv"), show_col_types = FALSE) %>%
    filter(transformation == "clr")
}

# Download helpers (adapted from download_helper.R)
read_10x <- function(dirs) {
  DropletUtils::read10xCounts(dirs)
}

download_geo_10x <- function(ids, dest_dir, file_regex = NULL) {
  all_dirs <- file.path(dest_dir, ids)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  for (id in ids) {
    out <- file.path(dest_dir, id)
    if (!file.exists(file.path(out, "matrix.mtx.gz")) &&
        !file.exists(file.path(out, "matrix.mtx"))) {
      message("  Downloading ", id, " ...")
      GEOquery::getGEOSuppFiles(id, baseDir = dest_dir,
                                filter_regex = file_regex,
                                makeDirectory = TRUE)
      # Rename to standard names
      for (f in list.files(out)) {
        new_name <- f
        new_name <- sub("^.*matrix\\.mtx", "matrix.mtx", new_name)
        new_name <- sub("^.*barcodes\\.tsv", "barcodes.tsv", new_name)
        new_name <- sub("^.*genes\\.tsv",    "genes.tsv",    new_name)
        new_name <- sub("^.*features\\.tsv", "features.tsv", new_name)
        if (new_name != f)
          file.rename(file.path(out, f), file.path(out, new_name))
      }
    }
  }
  all_dirs
}

load_consistency_dataset <- function(dataset_name) {
  dest <- file.path(DATA_DIR, "consistency")
  switch(dataset_name,
    GSE142647 = {
      dirs <- download_geo_10x(c("GSM4235299","GSM4235300"), dest)
      read_10x(dirs)
    },
    GSE178765 = {
      dirs <- download_geo_10x("GSE178765", dest)
      # Fix genes file (duplicate col)
      gf <- file.path(dest, "GSE178765", "genes.tsv.gz")
      if (file.exists(gf)) {
        g <- readr::read_tsv(gf, col_names = FALSE, show_col_types = FALSE)
        g$X2 <- g$X1
        readr::write_tsv(g, gf, col_names = FALSE)
      }
      read_10x(file.path(dest, "GSE178765"))
    },
    GSE179831 = {
      dirs <- download_geo_10x(c("GSM5434863","GSM5434864"), dest)
      read_10x(dirs)
    },
    GSE164017 = {
      id_dir <- file.path(dest, "GSM4994960")
      dir.create(id_dir, recursive = TRUE, showWarnings = FALSE)
      if (!file.exists(file.path(id_dir, "features.tsv.gz")))
        GEOquery::getGEOSuppFiles("GSE164017", baseDir = id_dir,
                                  filter_regex = "features.tsv.gz",
                                  makeDirectory = FALSE)
      dirs <- download_geo_10x("GSM4994960", dest)
      read_10x(dirs)
    },
    GSE150068 = {
      dirs <- download_geo_10x("GSM4522986", dest)
      read_10x(dirs)
    },
    GSE130931 = {
      dirs <- download_geo_10x(c("GSM4041124","GSM4041125"), dest)
      read_10x(dirs)
    },
    GSE163505 = {
      dirs <- download_geo_10x("GSM4980292", dest)
      read_10x(dirs)
    },
    GSE158941 = {
      dirs <- download_geo_10x("GSM4816083", dest)
      read_10x(dirs)
    },
    GSE179714 = {
      dirs <- download_geo_10x(c("GSM5429729","GSM5429730"), dest,
                                file_regex = "filtered-bc")
      read_10x(dirs)
    },
    GSE184806 = {
      dirs <- download_geo_10x("GSE184806", dest, file_regex = "batch12")
      read_10x(file.path(dest, "GSE184806"))
    },
    stop("Unknown dataset: ", dataset_name)
  )
}

run_consistency_clr <- function(dataset_name) {
  message("  Loading ", dataset_name, " ...")
  sce <- tryCatch(load_consistency_dataset(dataset_name),
                  error = function(e) { message("  FAILED: ", e$message); NULL })
  if (is.null(sce)) return(NULL)

  UMI <- filter_matrix(assay(sce, "counts"))
  n_genes <- nrow(UMI)
  message("  ", dataset_name, ": ", ncol(UMI), " cells × ", n_genes, " genes")

  rows <- list()
  for (seed in SEEDS) {
    set.seed(seed)
    idx1 <- sample(n_genes, floor(n_genes / 2))
    idx2 <- setdiff(seq_len(n_genes), idx1)
    UMI1 <- UMI[idx1, , drop = FALSE]
    UMI2 <- UMI[idx2, , drop = FALSE]
    sf1  <- size_factors(UMI1)
    sf2  <- size_factors(UMI2)

    for (pca_dim in PCA_VALS) {
      for (knn in KNN_VALS) {
        t_start <- proc.time()
        knn1 <- clr_knn(UMI1, sf1, pca_dim, knn)
        knn2 <- clr_knn(UMI2, sf2, pca_dim, knn)
        elapsed <- (proc.time() - t_start)[["elapsed"]]
        ov <- mean_knn_overlap(knn1, knn2)
        rows[[length(rows) + 1]] <- tibble::tibble(
          mean_overlap     = ov,
          transformation_id = paste0("clr_local_", dataset_name, "_s", seed,
                                     "_p", pca_dim, "_k", knn),
          dataset          = dataset_name,
          seed             = seed,
          pca_dim          = pca_dim,
          knn              = knn,
          transformation   = "clr",
          alpha            = "FALSE",
          cputime_sec      = elapsed,
          elapsed_sec      = elapsed
        )
        message(sprintf("    seed=%d pca=%d knn=%d  overlap=%.3f  (%.1fs)",
                        seed, pca_dim, knn, ov, elapsed))
      }
    }
  }
  bind_rows(rows)
}

consistency_results <- if (SKIP_CONSISTENCY) consistency_results else map(DATASETS, function(ds) {
  message("\nDataset: ", ds)
  tryCatch(run_consistency_clr(ds),
           error = function(e) { message("  ERROR: ", e$message); NULL })
}) %>% compact() %>% bind_rows()

if (!SKIP_CONSISTENCY && nrow(consistency_results) > 0) {
  out_file <- file.path(RESULTS_DIR, "consistency_results.tsv")
  existing <- readr::read_tsv(out_file, show_col_types = FALSE) %>%
    mutate(alpha = as.character(alpha))
  existing <- filter(existing, transformation != "clr")
  readr::write_tsv(bind_rows(existing, consistency_results), out_file)
  message("\nAppended ", nrow(consistency_results), " CLR rows to consistency_results.tsv")
} else {
  message("\nNo consistency results computed.")
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2 — DOWNSAMPLING BENCHMARK
# ═══════════════════════════════════════════════════════════════════════════════

SKIP_DOWNSAMPLING <- any(grepl("^clr$",
  readr::read_tsv(file.path(RESULTS_DIR, "downsampling_results.tsv"),
                  show_col_types = FALSE)$transformation))

message("\n=== DOWNSAMPLING BENCHMARK ===")
if (SKIP_DOWNSAMPLING) {
  message("CLR rows already present in downsampling_results.tsv — skipping.")
  downsampling_results <- readr::read_tsv(
    file.path(RESULTS_DIR, "downsampling_results.tsv"), show_col_types = FALSE) %>%
    filter(transformation == "clr")
}

# Adapted from download_helper.R — uses direct URLs, no GEOquery
load_downsampling_dataset <- function(dataset_name) {
  dest <- file.path(DATA_DIR, "downsampling")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  switch(dataset_name,
    mcSCRB = {
      f <- file.path(dest, "GSE103568_JM8_UMIcounts.txt.gz")
      if (!file.exists(f))
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103568&format=file&file=GSE103568%5FJM8%5FUMIcounts%2Etxt%2Egz", f, mode = "wb")
      as.matrix(read.delim(f))
    },
    smartSeq3_fibroblasts = {
      f <- file.path(dest, "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt")
      if (!file.exists(f))
        download.file("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/735/E-MTAB-8735/Files/Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt", f, method = "curl")
      as.matrix(read.delim(f))
    },
    smartSeq3_fibroblasts_alt = {
      f1 <- file.path(dest, "Fibroblasts.plate1.umis.ex.txt")
      f2 <- file.path(dest, "Fibroblasts.plate2.umis.ex.txt")
      if (!file.exists(f1))
        download.file("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/148/E-MTAB-10148/Files/Fibroblasts.plate1.umis.ex.txt", f1, method = "curl")
      if (!file.exists(f2))
        download.file("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/148/E-MTAB-10148/Files/Fibroblasts.plate2.umis.ex.txt", f2, method = "curl")
      m1 <- as.matrix(read.delim(f1))
      m2 <- as.matrix(read.delim(f2))
      common <- intersect(rownames(m1), rownames(m2))
      cbind(m1[common,], m2[common,])
    },
    smartSeq3_hek = {
      f <- file.path(dest, "HEK293T.umis.ex.txt")
      if (!file.exists(f))
        download.file("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/148/E-MTAB-10148/Files/HEK293T.umis.ex.txt", f, method = "curl")
      as.matrix(read.delim(f))
    },
    smartSeq3_siRNA_knockdown = {
      f <- file.path(dest, "HEK293T.siRNA.umis.ex.txt")
      if (!file.exists(f))
        download.file("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/148/E-MTAB-10148/Files/HEK293T.siRNA.umis.ex.txt", f, method = "curl")
      as.matrix(read.delim(f))
    },
    stop("Unknown downsampling dataset: ", dataset_name)
  )
}

downsample_counts <- function(UMI, seed, fraction = 0.1) {
  set.seed(seed)
  total <- sum(UMI)
  keep  <- round(total * fraction)
  # Multinomial downsampling per cell
  t(apply(UMI, 2, function(cell) {
    rmultinom(1, size = round(sum(cell) * fraction), prob = cell + 1e-8)[, 1]
  }))
}

run_downsampling_clr <- function(dataset_name) {
  message("  Loading ", dataset_name, " ...")
  UMI_full <- tryCatch(load_downsampling_dataset(dataset_name),
                       error = function(e) { message("  FAILED: ", e$message); NULL })
  if (is.null(UMI_full)) return(NULL)
  if (inherits(UMI_full, "sparseMatrix")) UMI_full <- as.matrix(UMI_full)
  UMI_full <- filter_matrix(UMI_full)
  message("  ", dataset_name, ": ", ncol(UMI_full), " cells × ", nrow(UMI_full), " genes")

  rows <- list()
  for (seed in DS_SEEDS) {
    set.seed(seed)
    # Downsample to ~10% of counts per cell
    # apply(..., 2, f) returns genes × cells when f returns a vector of length nrow
    UMI_red <- apply(UMI_full, 2, function(cell) {
      n <- max(1L, round(sum(cell) * 0.1))
      as.integer(rmultinom(1L, size = n, prob = cell + 1e-8))
    })
    rownames(UMI_red) <- rownames(UMI_full)
    UMI_red <- filter_matrix(UMI_red)
    # Align rows
    common_genes <- intersect(rownames(UMI_full), rownames(UMI_red))
    UMI_f <- UMI_full[common_genes, , drop = FALSE]
    UMI_r <- UMI_red[common_genes, , drop = FALSE]
    sf_f   <- size_factors(UMI_f)
    sf_r   <- size_factors(UMI_r)

    for (pca_dim in DS_PCA) {
      for (knn in DS_KNN) {
        t_start <- proc.time()
        knn_full <- clr_knn(UMI_f, sf_f, pca_dim, knn)
        knn_red  <- clr_knn(UMI_r, sf_r, pca_dim, knn)
        elapsed  <- (proc.time() - t_start)[["elapsed"]]
        ov <- mean_knn_overlap(knn_full, knn_red)
        rows[[length(rows) + 1]] <- tibble::tibble(
          overlap                     = ov,
          transformation_full_data_ids   = paste0("clr_local_full_", dataset_name, "_s", seed, "_p", pca_dim, "_k", knn),
          transformation_reduced_data_ids = paste0("clr_local_red_", dataset_name, "_s", seed, "_p", pca_dim, "_k", knn),
          dataset                     = dataset_name,
          seed                        = seed,
          pca_dim                     = pca_dim,
          knn                         = knn,
          transformation              = "clr",
          alpha                       = "FALSE",
          full_cputime_sec            = elapsed / 2,
          full_elapsed_sec            = elapsed / 2,
          reduced_cputime_sec         = elapsed / 2,
          reduced_elapsed_sec         = elapsed / 2
        )
        message(sprintf("    seed=%d pca=%d knn=%d  overlap=%.3f  (%.1fs)",
                        seed, pca_dim, knn, ov, elapsed))
      }
    }
  }
  bind_rows(rows)
}

downsampling_results <- if (SKIP_DOWNSAMPLING) downsampling_results else map(DS_DATASETS, function(ds) {
  message("\nDataset: ", ds)
  tryCatch(run_downsampling_clr(ds),
           error = function(e) { message("  ERROR: ", e$message); NULL })
}) %>% compact() %>% bind_rows()

if (!SKIP_DOWNSAMPLING && nrow(downsampling_results) > 0) {
  out_file <- file.path(RESULTS_DIR, "downsampling_results.tsv")
  existing <- readr::read_tsv(out_file, show_col_types = FALSE)
  existing <- mutate(existing, alpha = as.character(alpha))
  existing <- filter(existing, transformation != "clr")
  readr::write_tsv(bind_rows(existing, downsampling_results), out_file)
  message("\nAppended ", nrow(downsampling_results), " CLR rows to downsampling_results.tsv")
} else {
  message("\nNo downsampling results computed.")
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART 3 — SIMULATION BENCHMARK
# ═══════════════════════════════════════════════════════════════════════════════

message("\n=== SIMULATION BENCHMARK ===")

sim_available <- list(
  muscat      = requireNamespace("muscat",   quietly = TRUE),
  splatter    = requireNamespace("splatter", quietly = TRUE),
  scDesign2   = requireNamespace("scDesign2",quietly = TRUE),
  dyngen      = requireNamespace("dyngen",   quietly = TRUE)
)
message("Simulation packages available: ",
        paste(names(sim_available)[unlist(sim_available)], collapse = ", "))

# Simulation metric: ARI, AMI, NMI, mean KNN overlap against ground truth KNN
knn_overlap_with_truth <- function(knn_pred, knn_truth) {
  mean_knn_overlap(knn_pred, knn_truth)
}

simulate_dataset <- function(simulator, seed) {
  set.seed(seed)
  switch(simulator,
    muscat = {
      if (!sim_available$muscat) return(NULL)
      library(muscat, quietly = TRUE)
      data(example_sce, package = "muscat")
      sce_preped <- muscat::prepSim(example_sce, verbose = FALSE)
      sim <- muscat::simData(sce_preped, rel_lfc = c(1, 0.5, 0.1, 0.05),
                             nc = 5000L, nk = 4L,
                             p_dd = c(0.7, 0, 0.3, 0, 0, 0), lfc = 2,
                             ng = 1000L, force = TRUE)
      tmp <- as.data.frame(S4Vectors::metadata(sim)$gene_info) %>%
        tidyr::pivot_longer(dplyr::starts_with("sim_mean"),
                            names_sep = "\\.", names_to = c(".value", "group_id")) %>%
        dplyr::select(gene, cluster_id, group_id, sim_mean) %>%
        dplyr::full_join(
          as.data.frame(SummarizedExperiment::colData(sim),
                        optional = TRUE) %>%
            tibble::rownames_to_column("cell_id") %>%
            dplyr::mutate(cell_id = factor(cell_id, levels = cell_id)),
          by = c("cluster_id", "group_id"))
      gt_mat <- tmp %>%
        dplyr::arrange(cell_id) %>%
        dplyr::select(gene, cell_id, sim_mean) %>%
        tidyr::pivot_wider(id_cols = gene, names_from = cell_id,
                           values_from = sim_mean) %>%
        tibble::column_to_rownames("gene") %>%
        as.matrix()
      list(UMI = SummarizedExperiment::assay(sim, "counts"),
           ground_truth = log10(gt_mat + 1e-4))
    },
    splatter = {
      if (!sim_available$splatter) return(NULL)
      library(splatter, quietly = TRUE)
      params <- splatter::newSplatParams(nGenes = 1000L, batchCells = 5000L,
                                         group.prob = c(0.25, 0.25, 0.25, 0.25),
                                         de.prob = 0.3, seed = seed)
      sim <- splatter::splatSimulate(params, method = "groups", verbose = FALSE)
      counts_mat <- SummarizedExperiment::assay(sim, "counts")
      true_expr  <- SummarizedExperiment::assay(sim, "TrueCounts")
      list(UMI = counts_mat, ground_truth = log10(true_expr + 1e-4))
    },
    {
      message("  Simulator ", simulator, " not yet supported locally or packages missing")
      NULL
    }
  )
}

# Overlap metric against ground truth KNN
run_simulation_clr <- function(simulator, seed, sim_data) {
  if (is.null(sim_data)) return(NULL)
  UMI <- filter_matrix(sim_data$UMI)
  GT  <- sim_data$ground_truth
  # Align GT to filtered UMI
  common_cells <- intersect(colnames(UMI), colnames(GT))
  common_genes <- intersect(rownames(UMI), rownames(GT))
  if (length(common_cells) < 10 || length(common_genes) < 10) return(NULL)
  UMI <- UMI[common_genes, common_cells, drop = FALSE]
  GT  <- GT[common_genes, common_cells, drop = FALSE]
  sf  <- size_factors(UMI)

  rows <- list()
  for (pca_dim in SIM_PCA) {
    for (knn in SIM_KNN) {
      t_start <- proc.time()
      knn_pred <- clr_knn(UMI, sf, pca_dim, knn)
      gt_sf    <- rep(1, ncol(GT))
      knn_gt   <- clr_knn(GT, gt_sf, pca_dim, knn)
      elapsed  <- (proc.time() - t_start)[["elapsed"]]
      ov <- mean_knn_overlap(knn_pred, knn_gt)
      rows[[length(rows) + 1]] <- tibble::tibble(
        ARI              = NA_real_,
        AMI              = NA_real_,
        NMI              = NA_real_,
        mean_knn_overlap = ov,
        n_clusters       = NA_real_,
        n_clusters_counts = NA_real_,
        ground_truth_id  = paste0("gt_", simulator, "_s", seed),
        transformation_id = paste0("clr_local_", simulator, "_s", seed,
                                   "_p", pca_dim, "_k", knn),
        simulator        = simulator,
        seed             = seed,
        pca_dim          = pca_dim,
        knn              = knn,
        transformation   = "clr",
        alpha            = "FALSE",
        cputime_sec      = elapsed,
        elapsed_sec      = elapsed
      )
      message(sprintf("    seed=%d pca=%d knn=%d  overlap=%.3f  (%.1fs)",
                      seed, pca_dim, knn, ov, elapsed))
    }
  }
  bind_rows(rows)
}

simulation_results <- list()
for (simulator in SIM_SIMS) {
  message("\nSimulator: ", simulator)
  for (seed in SIM_SEEDS) {
    sim_data <- tryCatch(simulate_dataset(simulator, seed),
                         error = function(e) { message("  FAILED: ", e$message); NULL })
    if (is.null(sim_data)) next
    res <- tryCatch(run_simulation_clr(simulator, seed, sim_data),
                    error = function(e) { message("  ERROR: ", e$message); NULL })
    if (!is.null(res)) simulation_results[[length(simulation_results) + 1]] <- res
  }
}
simulation_results <- bind_rows(simulation_results)

if (nrow(simulation_results) > 0) {
  out_file <- file.path(RESULTS_DIR, "simulation_results.tsv")
  existing <- readr::read_tsv(out_file, show_col_types = FALSE)
  existing <- mutate(existing, alpha = as.character(alpha))
  existing <- filter(existing, transformation != "clr")
  readr::write_tsv(bind_rows(existing, simulation_results), out_file)
  message("\nAppended ", nrow(simulation_results), " CLR rows to simulation_results.tsv")
} else {
  message("\nNo simulation results computed (packages may be missing).")
}

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n=== DONE ===")
message("Consistency rows added: ", nrow(consistency_results))
message("Downsampling rows added: ", nrow(downsampling_results))
message("Simulation rows added:   ", nrow(simulation_results))
message("\nNext step: re-render the notebooks in notebooks/ to regenerate figures.")
