library(tidyverse)
library(SingleCellExperiment)
library(MatrixGenerics)
setwd("/Users/lpachter/Dropbox/claude/projects/edgetemp/transformGamPoi-Paper/notebooks")
source("annotation_helper.R")

# ── 1. Download Svensson 2017 data (pure technical noise, empty droplets) ─────
h5ad_path <- "../extra_data/svensson_2017_1.h5ad"
if (!file.exists(h5ad_path)) {
  cat("Downloading Svensson 2017 dataset...\n")
  download.file(
    "https://data.caltech.edu/records/6fy1w-dgm81/files/Svensson%20et%20al%202017%20(1).h5ad?download=1",
    h5ad_path
  )
}

set.seed(1)
se <- zellkonverter::readH5AD(h5ad_path, reader = "R")
se <- se[, sample.int(ncol(se), 400)]
cat("Dimensions:", nrow(se), "genes x", ncol(se), "cells\n")

size_factors <- Matrix::colSums(assay(se, "X"))
size_factors <- size_factors / mean(size_factors)
Y <- as.matrix(assay(se, "X"))
Y <- Y[Matrix::rowSums(Y) > 0, ]
cat("After filtering zeros:", nrow(Y), "genes\n")

# ── 2. Define transformations ─────────────────────────────────────────────────
clr_transform <- function(UMI, sf) {
  logpf <- log1p(sweep(UMI, 2L, sf, "/"))
  sweep(logpf, 2L, colMeans(logpf), "-")
}
clr_alpha_transform <- function(UMI, sf, alpha = 0.05) {
  logpf <- log(sweep(UMI, 2L, sf, "/") + 1 / (4 * alpha))
  sweep(logpf, 2L, colMeans(logpf), "-")
}

logp1_transform   <- function(UMI, sf) log1p(sweep(UMI, 2L, sf, "/"))
logp_cpm          <- function(UMI, sf) log1p(sweep(UMI, 2L, sf / 1e6, "/"))
logp1_size_normed <- function(UMI, sf) {
  lp <- log1p(sweep(UMI, 2L, sf, "/"))
  sweep(lp, 2L, colMeans(lp) / mean(colMeans(lp)), "/")
}
acosh_transform <- function(UMI, sf, alpha = 0.05) {
  # acosh(2*alpha*y/s + 1) / (2*sqrt(alpha))  [delta method variance stabilization]
  acosh(2 * alpha * sweep(UMI, 2L, sf, "/") + 1) / (2 * sqrt(alpha))
}
logp_alpha_transform <- function(UMI, sf, alpha = 0.05) {
  transformGamPoi::shifted_log_transform(UMI, size_factors = sf,
                                         pseudo_count = 1 / (4 * alpha))
}
pearson_clip_transform <- function(UMI, sf) {
  # Implement clipped Pearson residuals directly (avoids matrixStats useNames bug)
  mu <- outer(rowMeans(UMI / matrix(sf, nrow=nrow(UMI), ncol=ncol(UMI), byrow=TRUE)),
              sf)  # genes x cells expected counts
  r  <- (UMI - mu) / sqrt(mu + mu^2 * 0.05)
  pmax(pmin(r, sqrt(ncol(UMI))), -sqrt(ncol(UMI)))  # clip at sqrt(n)
}
raw_counts        <- function(UMI, sf) UMI
scaled_raw_counts <- function(UMI, sf) sweep(UMI, 2L, sf, "/")

# ── 3. Apply transformations ──────────────────────────────────────────────────
cat("Applying transformations...\n")
transformed_dat <- list()

safe_apply <- function(name, fn) {
  cat(" -", name, "... ")
  result <- tryCatch(fn(), error = function(e) { cat("SKIPPED:", conditionMessage(e), "\n"); NULL })
  if (!is.null(result)) cat("done\n")
  result
}

transformed_dat[["clr"]]               <- safe_apply("clr",               function() clr_transform(Y, size_factors))
transformed_dat[["clr_alpha"]]         <- safe_apply("clr_alpha",         function() clr_alpha_transform(Y, size_factors))
transformed_dat[["logp1"]]             <- safe_apply("logp1",             function() logp1_transform(Y, size_factors))
transformed_dat[["logp_cpm"]]          <- safe_apply("logp_cpm",          function() logp_cpm(Y, size_factors))
transformed_dat[["logp1_size_normed"]] <- safe_apply("logp1_size_normed", function() logp1_size_normed(Y, size_factors))
transformed_dat[["acosh"]]             <- safe_apply("acosh",             function() acosh_transform(Y, size_factors))
transformed_dat[["logp_alpha"]]        <- safe_apply("logp_alpha",        function() logp_alpha_transform(Y, size_factors))
transformed_dat[["pearson_clip"]]      <- safe_apply("pearson_clip",      function() pearson_clip_transform(Y, size_factors))
transformed_dat[["raw_counts"]]        <- safe_apply("raw_counts",        function() raw_counts(Y, size_factors))
transformed_dat[["scaled_raw_counts"]] <- safe_apply("scaled_raw_counts", function() scaled_raw_counts(Y, size_factors))

# Remove NULLs
transformed_dat <- Filter(Negate(is.null), transformed_dat)
cat("Transformations applied:", paste(names(transformed_dat), collapse=", "), "\n")

# ── 4. PCA on each transformation ─────────────────────────────────────────────
cat("Running PCA...\n")
pca_results <- map(names(transformed_dat), function(name) {
  cat(" -", name, "...\n")
  dat <- transformed_dat[[name]]
  # Ensure finite values
  dat[!is.finite(dat)] <- 0
  irlba::prcomp_irlba(t(dat), n = 10)$x
}) %>% set_names(names(transformed_dat))

# ── 5. Canonical correlation with size factor ─────────────────────────────────
cat("Computing canonical correlations...\n")
can_cor_res <- map_dfr(names(pca_results), function(name) {
  tibble(
    transformation = name,
    canonical_correlation = cancor(pca_results[[name]], size_factors)$cor[1]
  )
}) %>%
  left_join(trans_families, by = "transformation") %>%
  arrange(canonical_correlation)

cat("\n=== Canonical correlation with size factor (lower = better) ===\n")
print(can_cor_res %>% select(transformation, family, canonical_correlation) %>% arrange(canonical_correlation), n = 20)

# ── 6. Plot ───────────────────────────────────────────────────────────────────
family_order <- c("negative_control", "count_model", "latent_expr", "glm_residual", "delta_method")
can_cor_res <- can_cor_res %>%
  mutate(
    family = factor(family, levels = family_order),
    label  = trans_labels_plain[transformation],
    transformation = factor(transformation, levels = can_cor_res$transformation[order(can_cor_res$canonical_correlation)])
  )

p <- ggplot(can_cor_res, aes(x = canonical_correlation, y = reorder(label, -canonical_correlation), fill = family)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = trans_families_colors, labels = trans_families_labels_long,
                    name = "Method family") +
  scale_x_continuous(expand = expansion(add = c(0, 0.05)), limits = c(0, 1)) +
  labs(
    x = "Canonical Correlation with Size Factor",
    y = NULL,
    title = "Size factor influence on PCA (Svensson 2017 empty droplets)",
    subtitle = "Lower = less depth confounding. Data is pure noise — any structure is artifactual."
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title.position = "plot", legend.position = "bottom")

ggsave("../output/size_factor_canonical_correlation_with_clr.pdf", p, width = 7, height = 5)
ggsave("../output/size_factor_canonical_correlation_with_clr.png", p, width = 7, height = 5, dpi = 150)
cat("\nSaved to output/size_factor_canonical_correlation_with_clr.pdf\n")

# Save numerical results
write_csv(can_cor_res %>% select(transformation, family, canonical_correlation),
          "../output/size_factor_canonical_correlation_with_clr.csv")
