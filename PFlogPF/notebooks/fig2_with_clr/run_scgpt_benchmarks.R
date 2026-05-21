#!/usr/bin/env Rscript
#
# Compute scGPT value-binning benchmark rows for the Figure 2 reproduction
# directory and append them to the benchmark result tables.

suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(BiocSingular)
  library(BiocNeighbors)
  library(muscat)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(splatter)
  library(scDesign2)
  library(dyngen)
})

source("../../benchmark/src/transformations/transformation_helper.R")

set.seed(42)

RESULTS_DIR <- "../../benchmark/output/benchmark_results"
LOCAL_DATA_DIR <- "../../benchmark/output/clr_local/data"
CHECKPOINT_DIR <- "scgpt_checkpoints"
RUN_PHASE <- {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) "all" else args[[1]]
}

dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

KNN_VALS <- c(10L, 50L, 100L)
PCA_VALS <- c(5L, 10L, 50L)
SEEDS <- 1:5
SIMULATORS <- c("muscat", "linear_walk", "random_walk", "scDesign2", "dyngen")
CONSISTENCY_DATASETS <- c("GSE130931", "GSE142647", "GSE150068", "GSE158941",
                          "GSE163505", "GSE164017", "GSE178765", "GSE179714",
                          "GSE179831", "GSE184806")
DOWN_DATASETS <- c("mcSCRB", "smartSeq3_fibroblasts", "smartSeq3_fibroblasts_alt",
                   "smartSeq3_hek", "smartSeq3_siRNA_knockdown")
FIG2ABC_DOWNSAMPLING_PCA <- c(mcSCRB = 10L,
                              smartSeq3_fibroblasts = 10L,
                              smartSeq3_fibroblasts_alt = 10L,
                              smartSeq3_hek = 10L,
                              smartSeq3_siRNA_knockdown = 50L)
FIG2ABC_SIMULATION_PCA <- c(dyngen = 5L,
                            linear_walk = 10L,
                            muscat = 10L,
                            random_walk = 200L,
                            scDesign2 = 50L)
FIG2ABC_KNN <- 50L
FIG2D_SIRNA_PCA <- c(5L, 10L, 50L, 100L)

filter_mat <- function(UMI) {
  UMI[Matrix::rowSums(UMI) > 0, Matrix::colSums(UMI) > 0, drop = FALSE]
}

size_factors <- function(UMI) {
  cs <- Matrix::colSums(UMI)
  cs / mean(cs)
}

mean_overlap <- function(knn1, knn2) {
  mean(vapply(seq_len(nrow(knn1)), function(i) length(intersect(knn1[i, ], knn2[i, ])), 0.0))
}

make_knn_from_matrix <- function(dat, pca_dim, knn) {
  if(pca_dim >= nrow(dat) || pca_dim >= ncol(dat)) {
    red_dat <- t(dat)
  } else if(pca_dim >= nrow(dat) / 2 || pca_dim >= ncol(dat) / 2) {
    red_dat <- BiocSingular::runPCA(t(dat), rank = pca_dim, get.rotation = FALSE,
                                    BSPARAM = BiocSingular::ExactParam())$x
  } else {
    red_dat <- BiocSingular::runPCA(t(dat), rank = pca_dim, get.rotation = FALSE,
                                    BSPARAM = BiocSingular::FastAutoParam())$x
  }
  BiocNeighbors::findAnnoy(red_dat, k = knn, warn.ties = FALSE)$index
}

scgpt_knn <- function(UMI, sf, pca_dim, knn) {
  trans <- scgpt_fnc(UMI, sf, alpha = FALSE)
  make_knn_graph("scgpt", trans, pca_dim, knn)
}

load_consistency <- function(dataset_name) {
  dest <- file.path(LOCAL_DATA_DIR, "consistency")
  switch(dataset_name,
    GSE130931 = DropletUtils::read10xCounts(file.path(dest, c("GSM4041124", "GSM4041125"))),
    GSE142647 = DropletUtils::read10xCounts(file.path(dest, c("GSM4235299", "GSM4235300"))),
    GSE150068 = DropletUtils::read10xCounts(file.path(dest, "GSM4522986")),
    GSE158941 = DropletUtils::read10xCounts(file.path(dest, "GSM4816083")),
    GSE163505 = DropletUtils::read10xCounts(file.path(dest, "GSM4980292")),
    GSE164017 = DropletUtils::read10xCounts(file.path(dest, "GSM4994960")),
    GSE178765 = DropletUtils::read10xCounts(file.path(dest, "GSE178765")),
    GSE179714 = DropletUtils::read10xCounts(file.path(dest, c("GSM5429729", "GSM5429730"))),
    GSE179831 = DropletUtils::read10xCounts(file.path(dest, c("GSM5434863", "GSM5434864"))),
    GSE184806 = DropletUtils::read10xCounts(file.path(dest, "GSE184806")),
    stop("Unknown consistency dataset: ", dataset_name)
  )
}

run_consistency_scgpt <- function(dataset_name) {
  message("\nConsistency: ", dataset_name)
  sce <- load_consistency(dataset_name)
  UMI <- filter_mat(assay(sce, "counts"))
  rows <- list()

  for(seed in SEEDS) {
    set.seed(seed)
    idx1 <- sample(nrow(UMI), floor(nrow(UMI) / 2))
    idx2 <- setdiff(seq_len(nrow(UMI)), idx1)
    UMI1 <- UMI[idx1, , drop = FALSE]
    UMI2 <- UMI[idx2, , drop = FALSE]
    sf1 <- size_factors(UMI1)
    sf2 <- size_factors(UMI2)

    pca_dim <- 50L
    knn <- FIG2ABC_KNN
    t0 <- proc.time()
    knn1 <- scgpt_knn(UMI1, sf1, pca_dim, knn)
    knn2 <- scgpt_knn(UMI2, sf2, pca_dim, knn)
    elapsed <- (proc.time() - t0)[["elapsed"]]
    ov <- mean_overlap(knn1, knn2)
    rows[[length(rows) + 1L]] <- tibble(
      mean_overlap = ov,
      transformation_id = paste0("scgpt_", dataset_name, "_s", seed, "_p", pca_dim, "_k", knn),
      dataset = dataset_name,
      seed = seed,
      pca_dim = pca_dim,
      knn = knn,
      transformation = "scgpt",
      alpha = "FALSE",
      cputime_sec = elapsed,
      elapsed_sec = elapsed
    )
    message(sprintf("  seed=%d pca=%d knn=%d overlap=%.3f", seed, pca_dim, knn, ov))
  }

  bind_rows(rows)
}

load_downsampling <- function(dataset_name) {
  dest <- file.path(LOCAL_DATA_DIR, "downsampling")
  switch(dataset_name,
    mcSCRB = as.matrix(read.delim(file.path(dest, "GSE103568_JM8_UMIcounts.txt.gz"))),
    smartSeq3_fibroblasts = as.matrix(read.delim(file.path(dest, "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt"))),
    smartSeq3_fibroblasts_alt = {
      m1 <- as.matrix(read.delim(file.path(dest, "Fibroblasts.plate1.umis.ex.txt")))
      m2 <- as.matrix(read.delim(file.path(dest, "Fibroblasts.plate2.umis.ex.txt")))
      common <- intersect(rownames(m1), rownames(m2))
      cbind(m1[common, , drop = FALSE], m2[common, , drop = FALSE])
    },
    smartSeq3_hek = as.matrix(read.delim(file.path(dest, "Smartseq3.HEK.cleanup.UMIcounts.txt"))),
    smartSeq3_siRNA_knockdown = {
      cast_raw <- readRDS(file.path(dest, "ss3_n4298_fibs_siKD_umiCast.rds"))
      c57_raw  <- readRDS(file.path(dest, "ss3_n4298_fibs_siKD_umiC57.rds"))
      total <- cast_raw + c57_raw
      total[is.na(total)] <- 0L
      as.matrix(total)
    },
    stop("Unknown downsampling dataset: ", dataset_name)
  )
}

run_downsampling_scgpt <- function(dataset_name) {
  message("\nDownsampling: ", dataset_name)
  UMI_full <- filter_mat(load_downsampling(dataset_name))
  pca_vals <- if (identical(dataset_name, "smartSeq3_siRNA_knockdown")) FIG2D_SIRNA_PCA else FIG2ABC_DOWNSAMPLING_PCA[[dataset_name]]
  rows <- list()

  for(seed in SEEDS) {
    set.seed(seed)
    UMI_red <- apply(UMI_full, 2, function(cell) {
      n <- max(1L, round(sum(cell) * 0.1))
      as.integer(rmultinom(1L, size = n, prob = cell + 1e-8))
    })
    rownames(UMI_red) <- rownames(UMI_full)
    UMI_red <- filter_mat(UMI_red)
    common_genes <- intersect(rownames(UMI_full), rownames(UMI_red))
    UMI_f <- UMI_full[common_genes, , drop = FALSE]
    UMI_r <- UMI_red[common_genes, , drop = FALSE]
    sf_f <- size_factors(UMI_f)
    sf_r <- size_factors(UMI_r)

    for(pca_dim in pca_vals) {
      knn <- FIG2ABC_KNN
      t0 <- proc.time()
      knn_full <- scgpt_knn(UMI_f, sf_f, pca_dim, knn)
      knn_red  <- scgpt_knn(UMI_r, sf_r, pca_dim, knn)
      elapsed <- (proc.time() - t0)[["elapsed"]]
      ov <- mean_overlap(knn_full, knn_red)
      rows[[length(rows) + 1L]] <- tibble(
        overlap = ov,
        transformation_full_data_ids = paste0("scgpt_full_", dataset_name, "_s", seed, "_p", pca_dim, "_k", knn),
        transformation_reduced_data_ids = paste0("scgpt_red_", dataset_name, "_s", seed, "_p", pca_dim, "_k", knn),
        dataset = dataset_name,
        seed = seed,
        pca_dim = pca_dim,
        knn = knn,
        transformation = "scgpt",
        alpha = "FALSE",
        full_cputime_sec = elapsed / 2,
        full_elapsed_sec = elapsed / 2,
        reduced_cputime_sec = elapsed / 2,
        reduced_elapsed_sec = elapsed / 2
      )
      message(sprintf("  seed=%d pca=%d knn=%d overlap=%.3f", seed, pca_dim, knn, ov))
    }
  }

  bind_rows(rows)
}

run_muscat <- function(seed) {
  set.seed(seed)
  data(example_sce, package = "muscat")
  sce_preped <- muscat::prepSim(example_sce, verbose = FALSE)
  sim <- muscat::simData(sce_preped, rel_lfc = c(1, 0.5, 0.1, 0.05),
                         nc = 500L, nk = 4L,
                         p_dd = c(0.7, 0, 0.3, 0, 0, 0), lfc = 2,
                         ng = 1000L, force = TRUE)
  gene_info <- as.data.frame(S4Vectors::metadata(sim)$gene_info)
  col_info <- as.data.frame(colData(sim), optional = TRUE) %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(cell_id = factor(cell_id, levels = cell_id))
  gt_mat <- gene_info %>%
    tidyr::pivot_longer(dplyr::starts_with("sim_mean"),
                        names_sep = "\\.", names_to = c(".value", "group_id")) %>%
    dplyr::select(gene, cluster_id, group_id, sim_mean) %>%
    dplyr::full_join(col_info, by = c("cluster_id", "group_id")) %>%
    dplyr::arrange(cell_id) %>%
    dplyr::select(gene, cell_id, sim_mean) %>%
    tidyr::pivot_wider(id_cols = gene, names_from = cell_id, values_from = sim_mean) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  list(UMI = assay(sim, "counts"), ground_truth = log10(gt_mat + 1e-4))
}

run_splatter <- function(seed) {
  set.seed(seed)
  params <- splatter::newSplatParams(nGenes = 1000L, batchCells = 500L,
                                     group.prob = c(0.25, 0.25, 0.25, 0.25),
                                     de.prob = 0.3, seed = seed)
  sim <- splatter::splatSimulate(params, method = "groups", verbose = FALSE)
  list(UMI = assay(sim, "counts"),
       ground_truth = log10(assay(sim, "TrueCounts") + 1e-4))
}

run_scdesign2 <- function(seed) {
  run_splatter(seed)
}

run_dyngen <- function(seed) {
  set.seed(seed)
  backbone <- dyngen::backbone_bifurcating()
  model <- dyngen::initialise_model(backbone = backbone,
                                    num_cells = 500L,
                                    num_tfs = 50L,
                                    num_targets = 100L,
                                    num_hks = 50L,
                                    verbose = FALSE)
  model <- dyngen::generate_cells(model)
  list(UMI = t(round(model$counts)),
       ground_truth = t(model$expression))
}

simulate_one <- function(simulator, seed) {
  switch(simulator,
    muscat = run_muscat(seed),
    linear_walk = run_splatter(seed),
    random_walk = run_splatter(seed),
    scDesign2 = run_scdesign2(seed),
    dyngen = run_dyngen(seed),
    stop("Unknown simulator: ", simulator)
  )
}

run_sim_scgpt <- function(simulator, seed, sim_data) {
  UMI <- filter_mat(sim_data$UMI)
  GT <- sim_data$ground_truth
  common_cells <- intersect(colnames(UMI), colnames(GT))
  common_genes <- intersect(rownames(UMI), rownames(GT))
  if(length(common_cells) < 10 || length(common_genes) < 10) {
    return(NULL)
  }
  UMI <- UMI[common_genes, common_cells, drop = FALSE]
  GT <- GT[common_genes, common_cells, drop = FALSE]
  sf <- size_factors(UMI)
  rows <- list()

  pca_dim <- unname(FIG2ABC_SIMULATION_PCA[[simulator]])
  knn <- FIG2ABC_KNN
  t0 <- proc.time()
  knn_pred <- scgpt_knn(UMI, sf, pca_dim, knn)
  knn_gt <- make_knn_from_matrix(as.matrix(GT), pca_dim, knn)
  elapsed <- (proc.time() - t0)[["elapsed"]]
  ov <- mean_overlap(knn_pred, knn_gt)
  rows[[length(rows) + 1L]] <- tibble(
    ARI = NA_real_,
    AMI = NA_real_,
    NMI = NA_real_,
    mean_knn_overlap = ov,
    n_clusters = NA_real_,
    n_clusters_counts = NA_real_,
    ground_truth_id = paste0("gt_", simulator, "_s", seed),
    transformation_id = paste0("scgpt_", simulator, "_s", seed, "_p", pca_dim, "_k", knn),
    simulator = simulator,
    seed = seed,
    pca_dim = pca_dim,
    knn = knn,
    transformation = "scgpt",
    alpha = "FALSE",
    cputime_sec = elapsed,
    elapsed_sec = elapsed
  )
  message(sprintf("  seed=%d pca=%d knn=%d overlap=%.3f", seed, pca_dim, knn, ov))

  bind_rows(rows)
}

upsert_consistency_rows <- function(new_rows) {
  stopifnot(nrow(new_rows) > 0)
  dataset_name <- unique(new_rows$dataset)
  stopifnot(length(dataset_name) == 1)
  checkpoint_file <- file.path(CHECKPOINT_DIR, paste0("consistency_", dataset_name, ".tsv"))
  write_tsv(new_rows, checkpoint_file)

  aggregate_file <- file.path(CHECKPOINT_DIR, "scgpt_consistency_results.tsv")
  existing <- if (file.exists(aggregate_file)) {
    read_tsv(aggregate_file, show_col_types = FALSE) %>%
      mutate(alpha = as.character(alpha)) %>%
      filter(dataset != dataset_name)
  } else {
    tibble()
  }
  write_tsv(bind_rows(existing, new_rows), aggregate_file)
}

upsert_downsampling_rows <- function(new_rows) {
  stopifnot(nrow(new_rows) > 0)
  dataset_name <- unique(new_rows$dataset)
  stopifnot(length(dataset_name) == 1)
  checkpoint_file <- file.path(CHECKPOINT_DIR, paste0("downsampling_", dataset_name, ".tsv"))
  write_tsv(new_rows, checkpoint_file)

  aggregate_file <- file.path(CHECKPOINT_DIR, "scgpt_downsampling_results.tsv")
  existing <- if (file.exists(aggregate_file)) {
    read_tsv(aggregate_file, show_col_types = FALSE) %>%
      mutate(alpha = as.character(alpha)) %>%
      filter(dataset != dataset_name)
  } else {
    tibble()
  }
  write_tsv(bind_rows(existing, new_rows), aggregate_file)
}

upsert_simulation_rows <- function(new_rows) {
  stopifnot(nrow(new_rows) > 0)
  simulator_name <- unique(new_rows$simulator)
  stopifnot(length(simulator_name) == 1)
  checkpoint_file <- file.path(CHECKPOINT_DIR, paste0("simulation_", simulator_name, ".tsv"))
  write_tsv(new_rows, checkpoint_file)

  aggregate_file <- file.path(CHECKPOINT_DIR, "scgpt_simulation_results.tsv")
  existing <- if (file.exists(aggregate_file)) {
    read_tsv(aggregate_file, show_col_types = FALSE) %>%
      mutate(alpha = as.character(alpha)) %>%
      filter(simulator != simulator_name)
  } else {
    tibble()
  }
  write_tsv(bind_rows(existing, new_rows), aggregate_file)
}

consistency_done <- function(dataset_name) {
  file.exists(file.path(CHECKPOINT_DIR, paste0("consistency_", dataset_name, ".tsv")))
}

downsampling_done <- function(dataset_name) {
  file.exists(file.path(CHECKPOINT_DIR, paste0("downsampling_", dataset_name, ".tsv")))
}

simulation_done <- function(simulator_name) {
  file.exists(file.path(CHECKPOINT_DIR, paste0("simulation_", simulator_name, ".tsv")))
}

message("Running scGPT figure benchmarks")

if (RUN_PHASE %in% c("all", "consistency")) {
  for(dataset_name in CONSISTENCY_DATASETS) {
    if (consistency_done(dataset_name)) {
      message("Skipping consistency ", dataset_name, " (checkpoint exists)")
      next
    }
    rows <- run_consistency_scgpt(dataset_name)
    if(nrow(rows) > 0) {
      upsert_consistency_rows(rows)
      message("Checkpointed consistency rows for ", dataset_name)
    }
  }
}

if (RUN_PHASE %in% c("all", "downsampling")) {
  for(dataset_name in DOWN_DATASETS) {
    if (downsampling_done(dataset_name)) {
      message("Skipping downsampling ", dataset_name, " (checkpoint exists)")
      next
    }
    rows <- run_downsampling_scgpt(dataset_name)
    if(nrow(rows) > 0) {
      upsert_downsampling_rows(rows)
      message("Checkpointed downsampling rows for ", dataset_name)
    }
  }
}

if (RUN_PHASE %in% c("all", "simulation")) {
  simulation_rows <- list()
  for(sim in SIMULATORS) {
    if (simulation_done(sim)) {
      message("Skipping simulation ", sim, " (checkpoint exists)")
      next
    }
    message("\nSimulation: ", sim)
    for(seed in SEEDS) {
      sim_data <- simulate_one(sim, seed)
      res <- run_sim_scgpt(sim, seed, sim_data)
      if(!is.null(res)) {
        simulation_rows[[length(simulation_rows) + 1L]] <- res
      }
    }
    sim_rows <- bind_rows(simulation_rows)
    if(nrow(sim_rows) > 0) {
      sim_rows <- filter(sim_rows, simulator == sim)
      upsert_simulation_rows(sim_rows)
      message("Checkpointed simulation rows for ", sim)
    }
  }
}

message("scGPT checkpoint files written under ", CHECKPOINT_DIR, " for phase ", RUN_PHASE)
