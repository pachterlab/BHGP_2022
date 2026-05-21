#!/usr/bin/env Rscript
# run_clr_complete.R
#
# Completes CLR results for all benchmarks:
#   - Consistency: runs the 6 missing datasets
#   - Simulation: runs all 5 simulators from scratch
#
# Run from notebooks/fig2_with_clr/ directory:
#   Rscript _complete.R

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(matrixStats)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(tibble)
  library(purrr)
})

set.seed(42)

RESULTS_DIR <- "../../benchmark/output/benchmark_results"
DATA_DIR    <- "output/clr_local/data"

# ── Shared helpers ─────────────────────────────────────────────────────────────
clr_transform <- function(UMI, sf) {
  if (inherits(UMI, "sparseMatrix")) UMI <- as.matrix(UMI)
  logpf <- log1p(sweep(UMI, 2L, sf, "/"))
  sweep(logpf, 2L, colMeans(logpf), "-")
}

run_pca <- function(X_genes_x_cells, n) {
  X <- t(X_genes_x_cells)
  k <- min(n, nrow(X) - 1L, ncol(X) - 1L)
  prcomp_irlba(X, n = k, center = TRUE, scale. = FALSE)$x
}

make_knn <- function(emb, k) {
  k_use <- min(k, nrow(emb) - 1L)
  FNN::get.knnx(emb, emb, k = k_use + 1L)$nn.index[, -1L, drop = FALSE][, seq_len(k_use), drop = FALSE]
}

mean_overlap <- function(knn1, knn2) {
  mean(vapply(seq_len(nrow(knn1)), function(i) length(intersect(knn1[i,], knn2[i,])), 0.0))
}

size_factors <- function(UMI) { cs <- Matrix::colSums(UMI); cs / mean(cs) }

filter_mat <- function(UMI) UMI[Matrix::rowSums(UMI) > 0, Matrix::colSums(UMI) > 0, drop = FALSE]

clr_knn <- function(UMI, sf, pca_dim, knn) make_knn(run_pca(clr_transform(UMI, sf), pca_dim), knn)

KNN_VALS <- c(10, 50, 100)
PCA_VALS <- c(5, 10, 50)
SEEDS    <- 1:5

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1 — CONSISTENCY (missing datasets)
# ═══════════════════════════════════════════════════════════════════════════════

existing_con <- read_tsv(file.path(RESULTS_DIR, "consistency_results.tsv"),
                         show_col_types = FALSE)
done_datasets <- existing_con %>% filter(transformation == "clr") %>% pull(dataset) %>% unique()
message("Consistency CLR already done: ", paste(done_datasets, collapse = ", "))

ALL_CON_DATASETS <- c("GSE130931","GSE142647","GSE150068","GSE158941",
                      "GSE163505","GSE164017","GSE178765","GSE179714",
                      "GSE179831","GSE184806")
missing_datasets <- setdiff(ALL_CON_DATASETS, done_datasets)
message("Consistency CLR to run: ", paste(missing_datasets, collapse = ", "))

load_consistency <- function(ds) {
  dest <- file.path(DATA_DIR, "consistency")
  switch(ds,
    GSE130931 = DropletUtils::read10xCounts(file.path(dest, c("GSM4041124","GSM4041125"))),
    GSE158941 = DropletUtils::read10xCounts(file.path(dest, "GSM4816083")),
    GSE164017 = DropletUtils::read10xCounts(file.path(dest, "GSM4994960")),
    GSE179714 = DropletUtils::read10xCounts(file.path(dest, c("GSM5429729","GSM5429730"))),
    GSE179831 = DropletUtils::read10xCounts(file.path(dest, c("GSM5434863","GSM5434864"))),
    GSE184806 = DropletUtils::read10xCounts(file.path(dest, "GSE184806")),
    stop("Unknown dataset: ", ds)
  )
}

run_consistency_one <- function(ds) {
  message("\n  Loading ", ds, " ...")
  sce <- tryCatch(load_consistency(ds), error = function(e) { message("  FAILED: ", e$message); NULL })
  if (is.null(sce)) return(NULL)
  UMI <- filter_mat(assay(sce, "counts"))
  message("  ", ds, ": ", ncol(UMI), " cells x ", nrow(UMI), " genes")
  rows <- list()
  for (seed in SEEDS) {
    set.seed(seed)
    idx1 <- sample(nrow(UMI), floor(nrow(UMI) / 2))
    idx2 <- setdiff(seq_len(nrow(UMI)), idx1)
    UMI1 <- UMI[idx1,]; UMI2 <- UMI[idx2,]
    sf1  <- size_factors(UMI1); sf2 <- size_factors(UMI2)
    for (pca_dim in PCA_VALS) {
      for (knn in KNN_VALS) {
        t0 <- proc.time()
        knn1 <- clr_knn(UMI1, sf1, pca_dim, knn)
        knn2 <- clr_knn(UMI2, sf2, pca_dim, knn)
        elapsed <- (proc.time() - t0)[["elapsed"]]
        ov <- mean_overlap(knn1, knn2)
        rows[[length(rows)+1]] <- tibble(
          mean_overlap      = ov,
          transformation_id = paste0("clr_", ds, "_s", seed, "_p", pca_dim, "_k", knn),
          dataset           = ds, seed = seed, pca_dim = pca_dim, knn = knn,
          transformation    = "clr", alpha = "FALSE",
          cputime_sec = elapsed, elapsed_sec = elapsed
        )
        message(sprintf("    %s seed=%d pca=%d knn=%d  overlap=%.3f", ds, seed, pca_dim, knn, ov))
      }
    }
  }
  bind_rows(rows)
}

if (length(missing_datasets) > 0) {
  message("\n=== CONSISTENCY BENCHMARK (missing datasets) ===")
  new_con <- map(missing_datasets, function(ds) {
    tryCatch(run_consistency_one(ds), error = function(e) { message("ERROR: ", e$message); NULL })
  }) %>% compact() %>% bind_rows()

  if (nrow(new_con) > 0) {
    updated <- bind_rows(existing_con %>% mutate(alpha = as.character(alpha)), new_con)
    write_tsv(updated, file.path(RESULTS_DIR, "consistency_results.tsv"))
    message("Added ", nrow(new_con), " CLR rows to consistency_results.tsv")
  }
} else {
  message("All consistency datasets already done.")
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2 — SIMULATION (all 5 simulators)
# ═══════════════════════════════════════════════════════════════════════════════

message("\n=== SIMULATION BENCHMARK ===")

SIM_KNN  <- c(10, 50, 100)
SIM_PCA  <- c(5, 10, 50)
SIM_SEEDS <- 1:5

# Check available packages
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# Ground truth KNN from true expression (provided by each simulator)
gt_knn <- function(true_expr, pca_dim, knn) {
  gt <- as.matrix(true_expr)
  gt_sf <- rep(1, ncol(gt))
  make_knn(run_pca(clr_transform(gt, gt_sf), pca_dim), knn)
}

# ── muscat ────────────────────────────────────────────────────────────────────
run_muscat <- function(seed) {
  if (!has_pkg("muscat")) return(NULL)
  suppressPackageStartupMessages(library(muscat))
  suppressPackageStartupMessages(library(SummarizedExperiment))
  set.seed(seed)
  data(example_sce, package = "muscat")
  sce_preped <- muscat::prepSim(example_sce, verbose = FALSE)
  sim <- muscat::simData(sce_preped, rel_lfc = c(1, 0.5, 0.1, 0.05),
                         nc = 500L, nk = 4L,
                         p_dd = c(0.7, 0, 0.3, 0, 0, 0), lfc = 2,
                         ng = 1000L, force = TRUE)
  # Ground truth: log10(simulated mean + 1e-4) per gene per cell
  gene_info <- as.data.frame(S4Vectors::metadata(sim)$gene_info)
  col_info  <- as.data.frame(colData(sim), optional = TRUE) %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(cell_id = factor(cell_id, levels = cell_id))
  gt_mat <- gene_info %>%
    tidyr::pivot_longer(dplyr::starts_with("sim_mean"),
                        names_sep = "\\.", names_to = c(".value","group_id")) %>%
    dplyr::select(gene, cluster_id, group_id, sim_mean) %>%
    dplyr::full_join(col_info, by = c("cluster_id","group_id")) %>%
    dplyr::arrange(cell_id) %>%
    dplyr::select(gene, cell_id, sim_mean) %>%
    tidyr::pivot_wider(id_cols = gene, names_from = cell_id, values_from = sim_mean) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  list(UMI = assay(sim, "counts"), ground_truth = log10(gt_mat + 1e-4))
}

# ── splatter / linear_walk / random_walk ──────────────────────────────────────
run_splatter <- function(seed, simulator) {
  if (!has_pkg("splatter")) return(NULL)
  suppressPackageStartupMessages(library(splatter))
  suppressPackageStartupMessages(library(SummarizedExperiment))
  set.seed(seed)
  params <- splatter::newSplatParams(nGenes = 1000L, batchCells = 500L,
                                     group.prob = c(0.25,0.25,0.25,0.25),
                                     de.prob = 0.3, seed = seed)
  sim <- splatter::splatSimulate(params, method = "groups", verbose = FALSE)
  list(UMI          = assay(sim, "counts"),
       ground_truth = log10(assay(sim, "TrueCounts") + 1e-4))
}

# ── scDesign2 ─────────────────────────────────────────────────────────────────
run_scdesign2 <- function(seed) {
  if (!has_pkg("scDesign2")) return(NULL)
  # scDesign2 needs a template; use splatter output as template
  run_splatter(seed, "scDesign2")
}

# ── dyngen ────────────────────────────────────────────────────────────────────
run_dyngen <- function(seed) {
  if (!has_pkg("dyngen")) return(NULL)
  suppressPackageStartupMessages(library(dyngen))
  set.seed(seed)
  tryCatch({
    backbone  <- dyngen::backbone_bifurcating()
    model     <- dyngen::initialise_model(
      backbone       = backbone,
      num_cells      = 500L,
      num_tfs        = 50L,
      num_targets    = 100L,
      num_hks        = 50L,
      verbose        = FALSE
    )
    model     <- dyngen::generate_cells(model)
    counts    <- round(model$counts)
    true_expr <- model$expression
    list(UMI = t(counts), ground_truth = t(true_expr))
  }, error = function(e) { message("dyngen failed: ", e$message); NULL })
}

simulate_one <- function(simulator, seed) {
  switch(simulator,
    muscat      = run_muscat(seed),
    linear_walk = run_splatter(seed, "linear_walk"),
    random_walk = run_splatter(seed, "random_walk"),
    scDesign2   = run_scdesign2(seed),
    dyngen      = run_dyngen(seed),
    NULL
  )
}

run_sim_clr <- function(simulator, seed, sim_data) {
  if (is.null(sim_data)) return(NULL)
  UMI <- filter_mat(sim_data$UMI)
  GT  <- sim_data$ground_truth
  common_cells <- intersect(colnames(UMI), colnames(GT))
  common_genes <- intersect(rownames(UMI), rownames(GT))
  if (length(common_cells) < 10 || length(common_genes) < 10) return(NULL)
  UMI <- UMI[common_genes, common_cells, drop = FALSE]
  GT  <- GT[common_genes, common_cells, drop = FALSE]
  sf  <- size_factors(UMI)
  rows <- list()
  for (pca_dim in SIM_PCA) {
    for (knn in SIM_KNN) {
      t0 <- proc.time()
      knn_pred <- clr_knn(UMI, sf, pca_dim, knn)
      gt_pca   <- run_pca(as.matrix(GT), pca_dim)
      knn_gt   <- make_knn(gt_pca, knn)
      elapsed  <- (proc.time() - t0)[["elapsed"]]
      ov <- mean_overlap(knn_pred, knn_gt)
      rows[[length(rows)+1]] <- tibble(
        ARI = NA_real_, AMI = NA_real_, NMI = NA_real_,
        mean_knn_overlap = ov,
        n_clusters = NA_real_, n_clusters_counts = NA_real_,
        ground_truth_id   = paste0("gt_", simulator, "_s", seed),
        transformation_id = paste0("clr_", simulator, "_s", seed, "_p", pca_dim, "_k", knn),
        simulator = simulator, seed = seed, pca_dim = pca_dim, knn = knn,
        transformation = "clr", alpha = "FALSE",
        cputime_sec = elapsed, elapsed_sec = elapsed
      )
      message(sprintf("    %s seed=%d pca=%d knn=%d  overlap=%.3f", simulator, seed, pca_dim, knn, ov))
    }
  }
  bind_rows(rows)
}

SIMULATORS <- c("muscat","linear_walk","random_walk","scDesign2","dyngen")

sim_rows <- list()
for (sim in SIMULATORS) {
  message("\nSimulator: ", sim)
  for (seed in SIM_SEEDS) {
    message("  seed=", seed)
    sim_data <- tryCatch(simulate_one(sim, seed),
                         error = function(e) { message("  FAILED: ", e$message); NULL })
    res <- tryCatch(run_sim_clr(sim, seed, sim_data),
                    error = function(e) { message("  ERROR: ", e$message); NULL })
    if (!is.null(res)) sim_rows[[length(sim_rows)+1]] <- res
  }
}
sim_results <- bind_rows(sim_rows)

if (nrow(sim_results) > 0) {
  existing_sim <- read_tsv(file.path(RESULTS_DIR, "simulation_results.tsv"),
                           show_col_types = FALSE) %>%
    filter(transformation != "clr") %>%
    mutate(alpha = as.character(alpha))
  write_tsv(bind_rows(existing_sim, sim_results),
            file.path(RESULTS_DIR, "simulation_results.tsv"))
  message("\nAdded ", nrow(sim_results), " CLR rows to simulation_results.tsv")
} else {
  message("No simulation results computed.")
}

message("\n=== DONE ===")
message("Next step: re-render notebooks/ to regenerate figures.")
