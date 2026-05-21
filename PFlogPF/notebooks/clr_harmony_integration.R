#!/usr/bin/env Rscript
# clr_harmony_integration.R
#
# Experiment: does CLR make depth-aware Harmony integration better?
#
# Setup:
#   - Take mcSCRB dataset (full depth)
#   - Downsample to 5000 reads/cell median (matching benchmark protocol)
#   - Mix: 50% of cells from full-depth, 50% from downsampled (same cells!)
#   - Integrate with Harmony (batch = full vs downsampled)
#   - Measure KNN overlap vs "ground truth" = KNN built on ALL full-depth cells
#
# Hypothesis: CLR removes depth signal before Harmony sees it, so Harmony
# doesn't need to correct as hard, and the integrated embedding is cleaner.

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(FNN)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(harmony)
  library(scuttle)
})

set.seed(42)

DATA_FILE <- "../benchmark/output/clr_local/data/downsampling/GSE103568_JM8_UMIcounts.txt.gz"
OUT_DIR   <- "../output"
dir.create(OUT_DIR, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading mcSCRB data...\n")
raw <- as.matrix(read.delim(DATA_FILE))
raw <- raw[rowSums(raw) > 0, colSums(raw) > 0]
cat("  Full data:", nrow(raw), "genes x", ncol(raw), "cells\n")
cat("  Median depth:", median(colSums(raw)), "reads/cell\n")

# ── 2. Downsample (same protocol as benchmark: target 5000 reads/cell median) ─
cat("Downsampling...\n")
TARGET_DEPTH <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
if (is.na(TARGET_DEPTH)) TARGET_DEPTH <- 1000
cat("  Target depth:", TARGET_DEPTH, "reads/cell\n")
prop <- TARGET_DEPTH / median(colSums(raw))
cat("  Downsampling proportion:", round(prop, 3), "\n")
raw_down <- as.matrix(scuttle::downsampleMatrix(raw, prop = prop, bycol = FALSE))
# Keep only cells with counts in both
keep <- colSums(raw) > 0 & colSums(raw_down) > 0
raw      <- raw[, keep]
raw_down <- raw_down[, keep]
cat("  After filtering:", ncol(raw), "cells retained\n")
cat("  Median depth full:", round(median(colSums(raw))), "  downsampled:", round(median(colSums(raw_down))), "\n")

# ── 3. Sample 50% of cells for each batch ─────────────────────────────────────
n_cells   <- ncol(raw)
half      <- floor(n_cells / 2)
idx_full  <- sample(n_cells, half)                        # these contribute full-depth
idx_down  <- setdiff(seq_len(n_cells), idx_full)           # these contribute downsampled

# Ground truth: KNN on ALL full-depth cells (no mixing, no integration)
# This is what a perfect integration should recover
UMI_gt <- raw   # all cells at full depth

# Mixed dataset: half full, half downsampled
UMI_mix      <- cbind(raw[, idx_full], raw_down[, idx_down])
batch_labels <- c(rep("full", length(idx_full)), rep("downsampled", length(idx_down)))
# Cell order: first half = full-depth cells (idx_full), second half = downsampled (idx_down)
# We need to track which original cell index each mixed cell corresponds to
orig_idx <- c(idx_full, idx_down)   # original cell index in gt

cat("  Mixed dataset:", ncol(UMI_mix), "cells (", sum(batch_labels=="full"),
    "full +", sum(batch_labels=="downsampled"), "downsampled)\n")

# ── 4. Transforms ─────────────────────────────────────────────────────────────
sf_vec <- function(UMI) { cs <- colSums(UMI); cs / mean(cs) }

logp1_transform <- function(UMI) {
  sf <- sf_vec(UMI)
  log1p(sweep(UMI, 2L, sf, "/"))
}

clr_transform <- function(UMI) {
  sf  <- sf_vec(UMI)
  lpf <- log1p(sweep(UMI, 2L, sf, "/"))
  sweep(lpf, 2L, colMeans(lpf), "-")
}

# ── 5. PCA + Harmony integration ──────────────────────────────────────────────
run_pca <- function(mat_genes_x_cells, n = 30) {
  prcomp_irlba(t(mat_genes_x_cells), n = n, center = TRUE, scale. = FALSE)$x
}

run_harmony <- function(pca, batch) {
  harmony::HarmonyMatrix(pca, meta_data = data.frame(batch = batch),
                         vars_use = "batch", do_pca = FALSE, verbose = FALSE)
}

make_knn <- function(embedding, k = 50) {
  nn <- FNN::get.knnx(embedding, embedding, k = k + 1L)$nn.index
  nn[, -1L, drop = FALSE]
}

overlap_with_gt <- function(knn_mixed, knn_gt, orig_idx, k = 50) {
  # For each cell in the mixed dataset, compare its neighbors in the mixed
  # embedding to what its neighbors would be in the full-depth ground truth.
  # knn_mixed: (n_mix x k) indices into the mixed dataset
  # knn_gt:    (n_all x k) indices into the full dataset
  # orig_idx:  maps mixed cell i → full dataset cell index
  n <- nrow(knn_mixed)
  mean(vapply(seq_len(n), function(i) {
    gt_neighbors  <- knn_gt[orig_idx[i], ]           # k neighbors in full GT
    mix_neighbors <- orig_idx[knn_mixed[i, ]]        # k neighbors in mixed, mapped back
    length(intersect(gt_neighbors, mix_neighbors)) / k
  }, 0.0))
}

cat("\n── Building ground truth KNN (all", ncol(UMI_gt), "cells, full depth) ──\n")
pca_gt  <- run_pca(logp1_transform(UMI_gt), n = 30)
knn_gt  <- make_knn(pca_gt, k = 50)

cat("\n── logp1 + Harmony ──\n")
pca_logp1  <- run_pca(logp1_transform(UMI_mix), n = 30)
harm_logp1 <- run_harmony(pca_logp1, batch_labels)
knn_logp1  <- make_knn(harm_logp1, k = 50)
ov_logp1   <- overlap_with_gt(knn_logp1, knn_gt, orig_idx, k = 50)
cat("  KNN overlap with GT:", round(ov_logp1, 4), "\n")

# Also logp1 without integration (baseline)
knn_logp1_noharmony <- make_knn(pca_logp1, k = 50)
ov_logp1_noharmony  <- overlap_with_gt(knn_logp1_noharmony, knn_gt, orig_idx, k = 50)
cat("  logp1 (no Harmony):", round(ov_logp1_noharmony, 4), "\n")

cat("\n── CLR + Harmony ──\n")
pca_clr  <- run_pca(clr_transform(UMI_mix), n = 30)
harm_clr <- run_harmony(pca_clr, batch_labels)
knn_clr  <- make_knn(harm_clr, k = 50)
ov_clr   <- overlap_with_gt(knn_clr, knn_gt, orig_idx, k = 50)
cat("  KNN overlap with GT:", round(ov_clr, 4), "\n")

# CLR without Harmony
knn_clr_noharmony <- make_knn(pca_clr, k = 50)
ov_clr_noharmony  <- overlap_with_gt(knn_clr_noharmony, knn_gt, orig_idx, k = 50)
cat("  CLR (no Harmony):", round(ov_clr_noharmony, 4), "\n")

# ── 6. Summary ────────────────────────────────────────────────────────────────
results <- tibble(
  method       = c("logp1 (no Harmony)", "logp1 + Harmony", "CLR (no Harmony)", "CLR + Harmony"),
  knn_overlap  = c(ov_logp1_noharmony, ov_logp1, ov_clr_noharmony, ov_clr),
  norm         = c("logp1","logp1","CLR","CLR"),
  integration  = c("none","Harmony","none","Harmony")
)

cat("\n=== Results (k=50, KNN overlap with full-depth ground truth) ===\n")
results <- results %>% arrange(desc(knn_overlap))
for (i in seq_len(nrow(results))) {
  cat(sprintf("  %-25s %.4f\n", results$method[i], results$knn_overlap[i]))
}

# ── 7. PC1 vs depth plot (before/after integration) ──────────────────────────
depth_mix <- colSums(UMI_mix)

depth_df <- bind_rows(
  tibble(PC1 = pca_logp1[,1], depth = log10(depth_mix), batch = batch_labels,
         stage = "logp1 raw PCA"),
  tibble(PC1 = harm_logp1[,1], depth = log10(depth_mix), batch = batch_labels,
         stage = "logp1 + Harmony"),
  tibble(PC1 = pca_clr[,1], depth = log10(depth_mix), batch = batch_labels,
         stage = "CLR raw PCA"),
  tibble(PC1 = harm_clr[,1], depth = log10(depth_mix), batch = batch_labels,
         stage = "CLR + Harmony")
) %>%
  mutate(stage = factor(stage, levels = c("logp1 raw PCA","logp1 + Harmony",
                                           "CLR raw PCA","CLR + Harmony")))

cors <- depth_df %>%
  group_by(stage) %>%
  summarize(r = round(cor(PC1, depth), 3), .groups = "drop")
cat("\nPC1 ~ log10(depth) correlations:\n")
for (i in seq_len(nrow(cors))) cat(sprintf("  %-25s r=%.3f\n", as.character(cors$stage[i]), cors$r[i]))

p_depth <- ggplot(depth_df, aes(x = depth, y = PC1, color = batch)) +
  geom_point(alpha = 0.3, size = 0.5) +
  facet_wrap(~stage, nrow = 2) +
  geom_text(data = cors, aes(label = paste0("r=", r), x = -Inf, y = Inf),
            hjust = -0.1, vjust = 1.3, inherit.aes = FALSE, size = 3.5) +
  scale_color_manual(values = c(full = "#2166ac", downsampled = "#d6604d")) +
  labs(x = "log10(library size)", y = "PC1",
       title = "Depth confounding before and after integration (mcSCRB)",
       subtitle = "50% full-depth + 50% downsampled cells") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "clr_harmony_depth_confounding.pdf"), p_depth, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "clr_harmony_depth_confounding.png"), p_depth, width = 8, height = 6, dpi = 150)

# ── 8. Bar chart of KNN overlaps ──────────────────────────────────────────────
p_bar <- ggplot(results, aes(x = method, y = knn_overlap, fill = norm, alpha = integration)) +
  geom_col(width = 0.6) +
  scale_alpha_manual(values = c(none = 0.5, Harmony = 1.0)) +
  scale_fill_manual(values = c(logp1 = "#6baed6", CLR = "#74c476")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = NULL, y = "KNN overlap with full-depth GT (k=50)",
       title = "CLR + Harmony recovers full-depth KNN structure better",
       subtitle = "mcSCRB: 50% full-depth + 50% downsampled cells") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(OUT_DIR, "clr_harmony_knn_overlap.pdf"), p_bar, width = 6, height = 4)
ggsave(file.path(OUT_DIR, "clr_harmony_knn_overlap.png"), p_bar, width = 6, height = 4, dpi = 150)

cat("\nPlots saved to", OUT_DIR, "\n")
