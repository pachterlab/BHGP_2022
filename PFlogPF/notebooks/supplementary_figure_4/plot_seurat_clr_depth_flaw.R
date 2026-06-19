#!/usr/bin/env Rscript
# plot_seurat_clr_depth_flaw.R
#
# Same joint-embedding depth-flaw test as 07_, but repeated over many
# downsampling seeds and summarised as side-by-side violins per method.
#
# Metric (per cell): twin-recovery score = fraction of the OTHER depth group
# that is farther from the cell than its own twin.
#   1.0  -> the cell's own downsampled/full twin is its nearest cross-depth cell
#   0.5  -> twin is at a random rank (embedding driven by depth, not identity)
# Pooled over all cells and seeds; each seed's mean is overlaid as a point.
#
# Run with the R that has Seurat on its path:
#   Rscript plot_seurat_clr_depth_flaw.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(irlba)
  library(ggplot2)
  library(cowplot)
  library(reticulate)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
PFLOGPF_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
DATA_DIR <- file.path(PFLOGPF_ROOT, "benchmark/output/clr_local/data/downsampling")
OUT_DIR <- file.path(PFLOGPF_ROOT, "output")
DATASET  <- "smartSeq3_fibroblasts"
FRACTION <- 0.10
N_PC     <- 10
N_SEEDS  <- 20
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load once.
f <- file.path(DATA_DIR, "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt")
UMI <- as.matrix(read.delim(f))
UMI <- UMI[rowSums(UMI) > 0, colSums(UMI) > 0, drop = FALSE]
N <- ncol(UMI)
group <- rep(c("full", "downsampled"), each = N)   # depth label per joint column
message(sprintf("%s: %d cells x %d genes; median depth %.0f",
                DATASET, N, nrow(UMI), median(colSums(UMI))))

scclr_python <- Sys.getenv("SC_CLR_PYTHON", unset = "")
if (nzchar(scclr_python)) {
  use_python(scclr_python, required = TRUE)
}
scclr_path <- Sys.getenv("SC_CLR_PATH", unset = "")
if (!nzchar(scclr_path)) {
  sibling_scclr <- file.path(PFLOGPF_ROOT, "..", "..", "scclr", "python")
  if (dir.exists(sibling_scclr)) {
    scclr_path <- normalizePath(sibling_scclr, mustWork = TRUE)
  }
}
if (nzchar(scclr_path)) {
  py_run_string(sprintf("import sys; sys.path.insert(0, %s)", shQuote(scclr_path)))
}
scclr <- import("scclr")
np <- import("numpy")
message("Using scclr from: ", scclr$`__file__`)

size_factors <- function(M) { cs <- colSums(M); cs / mean(cs) }
estimate_alpha <- function(M) {
  sf <- size_factors(M)
  norm <- sweep(M, 2, sf, "/")
  mu <- rowMeans(norm)
  keep <- mu > 0
  norm <- norm[keep, , drop = FALSE]
  mu <- mu[keep]
  v <- apply(norm, 1, var)
  poisson_term <- mu * mean(1 / sf)
  alpha_g <- (v - poisson_term) / (mu^2)
  alpha_g <- alpha_g[is.finite(alpha_g) & alpha_g > 0]
  if (!length(alpha_g)) stop("could not estimate a positive overdispersion")
  median(alpha_g)
}
scclr_pflog <- function(M) {
  X <- np$array(t(M), dtype = "float64")
  sclr <- scclr$normalize(X, target = "mean", center = TRUE)
  t(sclr$to_dense())
}
transform_fns <- list(
  `raw counts`    = function(M) M,
  `Seurat "CLR"`  = function(M) as.matrix(NormalizeData(M, normalization.method = "CLR",
                                                        margin = 2, verbose = FALSE)),
  `log1pPF`       = function(M) log1p(sweep(M, 2, size_factors(M), "/")),
  `Ahlmann-Eltze–Huber` = function(M) {
    alpha <- estimate_alpha(M)
    log(sweep(M, 2, size_factors(M), "/") + 1 / (4 * alpha))
  },
  `PFlog` = scclr_pflog
)
ord <- names(transform_fns)

run_pca <- function(X_genes_x_cells, n) {
  X <- t(X_genes_x_cells)
  prcomp_irlba(X, n = min(n, nrow(X) - 1L, ncol(X) - 1L),
               center = TRUE, scale. = FALSE)$x
}

# Per-cell twin-recovery score from a joint embedding (first N rows = full,
# last N rows = downsampled twins, paired by index).
twin_recovery_scores <- function(emb, N) {
  Femb <- emb[seq_len(N), , drop = FALSE]
  Demb <- emb[N + seq_len(N), , drop = FALSE]
  # cross distance^2 between every full cell (row) and every down cell (col)
  D2 <- outer(rowSums(Femb^2), rowSums(Demb^2), "+") - 2 * Femb %*% t(Demb)
  diagd <- diag(D2)                                   # full i  <-> down i (twins)
  rank_full <- vapply(seq_len(N), function(i) sum(D2[i, ] <= diagd[i]), 0L)
  rank_down <- vapply(seq_len(N), function(j) sum(D2[, j] <= diagd[j]), 0L)
  1 - (c(rank_full, rank_down) - 1) / (N - 1)         # 1 = best, 0.5 = random
}

# Run over seeds.
rows <- list(); r2_rows <- list()
for (seed in seq_len(N_SEEDS)) {
  set.seed(seed)
  UMI_down <- apply(UMI, 2, function(cell) {
    n <- max(1L, round(sum(cell) * FRACTION))
    as.integer(rmultinom(1L, size = n, prob = cell + 1e-8))
  })
  rownames(UMI_down) <- rownames(UMI)
  joint <- cbind(UMI, UMI_down)
  joint <- joint[rowSums(joint) > 0, , drop = FALSE]

  for (m in ord) {
    set.seed(1000 + seed)
    emb <- run_pca(transform_fns[[m]](joint), N_PC)
    sc  <- twin_recovery_scores(emb, N)
    rows[[length(rows) + 1]] <- data.frame(method = m, seed = seed, score = sc)
    # fraction of leading-PC variance attributable to depth group (point-biserial R^2)
    r2_rows[[length(r2_rows) + 1]] <- data.frame(
      method = m, seed = seed,
      pc1_depth_R2 = summary(lm(emb[, 1] ~ group))$r.squared)
  }
  message(sprintf("seed %d/%d done", seed, N_SEEDS))
}
df <- do.call(rbind, rows)
df$method <- factor(df$method, levels = ord)
r2_df <- do.call(rbind, r2_rows)
r2_df$method <- factor(r2_df$method, levels = ord)

seed_means <- aggregate(score ~ method + seed, df, mean)
overall    <- aggregate(score ~ method, df, mean)
message("\nMean twin-recovery score by method (pooled over cells & seeds):")
print(overall[order(-overall$score), ], row.names = FALSE)
message("\nMean PC1-variance-from-depth (R^2) by method (over seeds):")
print(aggregate(pc1_depth_R2 ~ method, r2_df, mean), row.names = FALSE)

# Violin figure.
# Ahlmann-Eltze & Huber palette (ColorBrewer Set2), distinct hue per method.
# Our method PFlog is highlighted in red (#E41A1C), matching the
# manuscript main_benchmark_fig.pdf. (see notebooks/annotation_helper.R)
fills <- c("raw counts" = "#7d7d7d", `Seurat "CLR"` = "#000000",
           "log1pPF" = "#66C2A5", "Ahlmann-Eltze–Huber" = "#FC8D62",
           "PFlog" = "#E41A1C")

p <- ggplot(df, aes(method, score, fill = method)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
  geom_violin(scale = "width", width = 0.85, color = NA, alpha = 0.85,
              show.legend = FALSE) +
  geom_jitter(data = seed_means, aes(method, score), width = 0.07, height = 0,
              size = 1.4, color = "grey20", alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
               fill = "white", color = "black", show.legend = FALSE) +
  scale_fill_manual(values = fills) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "twin-recovery score") +
  theme_cowplot(font_size = 12) +
  guides(fill = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(OUT_DIR, "fig_seurat_clr_depth_flaw_violin.png"), p, width = 8, height = 5.5,
       dpi = 300, bg = "white")
message("\nSaved: ", file.path(OUT_DIR, "fig_seurat_clr_depth_flaw_violin.png"))
write.csv(seed_means, file.path(SCRIPT_DIR, "seurat_clr_depth_flaw_violin_seedmeans.csv"), row.names = FALSE)

# PC1-variance-from-depth violin: one R^2 per seed; lower = depth removed.
p_r2 <- ggplot(r2_df, aes(method, pc1_depth_R2, fill = method)) +
  geom_violin(scale = "width", width = 0.85, color = NA, alpha = 0.85,
              show.legend = FALSE) +
  geom_jitter(width = 0.07, height = 0, size = 1.4, color = "grey20", alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
               fill = "white", color = "black", show.legend = FALSE) +
  scale_fill_manual(values = fills) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = expression("PC1 variance from depth ("*R^2*")")) +
  theme_cowplot(font_size = 12) +
  guides(fill = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(OUT_DIR, "fig_seurat_clr_depth_pc1var_violin.png"), p_r2, width = 8, height = 5.5,
       dpi = 300, bg = "white")
message("Saved: ", file.path(OUT_DIR, "fig_seurat_clr_depth_pc1var_violin.png"))
write.csv(r2_df, file.path(SCRIPT_DIR, "seurat_clr_depth_pc1var_seeds.csv"), row.names = FALSE)

# Combined 2-panel figure: violin version of the bar figure.
# a: per-cell twin-recovery score   b: per-seed PC1-variance-from-depth (R^2)
fig2 <- plot_grid(p, p_r2, nrow = 1, labels = c("a", "b"))
ggsave(file.path(OUT_DIR, "seurat_clr_depth_flaw.png"), fig2, width = 9.5, height = 4.5,
       dpi = 300, bg = "white")
message("Saved: ", file.path(OUT_DIR, "seurat_clr_depth_flaw.png"))
