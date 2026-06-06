#!/usr/bin/env Rscript
# Figure 6: Pseudocount simulation.
#
# Coefficient of variation (CV) of the per-gene variances, plotted against the
# "correction factor" cf for a range of negative-binomial overdispersions alpha.
# Genes are transformed with the shifted log  log(y/sf + c),  pseudocount
# c = 1/(cf * alpha); cf = 4 is the delta-method / Anscombe optimum c = 1/(4*alpha).
# A well-chosen pseudocount stabilises the variance (uniform per-gene variances ->
# low CV); too large (small cf) leaves the raw variance spread, too small (large cf)
# over-amplifies low-count genes -> the CV rises at both ends.
#
# NOTE: this is a reproduction of the original Fig 6 (whose generating script was
# not in the repo). The simulation parameters below (gene-mean distribution,
# n_genes/n_cells, replicates) were tuned to match the published panels closely;
# they are not the exact originals, so absolute CV values are approximate (~1.5x)
# while the shapes and the optimum location match.
#
# Usage:  Rscript pseudocount_simulation.R   (writes ../output/pseudocount_simulation.{pdf,png})

set.seed(1)

n_genes <- 2000L
n_cells <- 500L
n_rep   <- 15L
alphas  <- c(0.05, 0.1, 0.3, 0.5, 1, 1.5)
cfs     <- c(2, 4, 6, 8, 12, 20)

# CV of per-gene variances after the shifted-log transform log(Y/sf + c).
cv_gene_var <- function(Y, sf, c) {
  Tr <- log(sweep(Y, 2L, sf, "/") + c)
  v  <- apply(Tr, 1L, var)            # per-gene variance across cells
  sd(v) / mean(v)
}

results <- vector("list", length(alphas)); names(results) <- as.character(alphas)
for (a in alphas) {
  M <- matrix(NA_real_, n_rep, length(cfs))
  for (r in seq_len(n_rep)) {
    gene_mean <- 10^runif(n_genes, 0.5, 1.1)              # gene means ~3..13
    s         <- rep(1, n_cells)                          # no size-factor variation
    mu        <- outer(gene_mean, s)
    Y         <- matrix(rnbinom(n_genes * n_cells, mu = mu, size = 1 / a),
                        n_genes, n_cells)                 # Var = mu + alpha*mu^2
    keep <- rowSums(Y) > 0
    Y    <- Y[keep, , drop = FALSE]
    sf   <- colSums(Y); sf <- sf / mean(sf)
    for (j in seq_along(cfs)) M[r, j] <- cv_gene_var(Y, sf, c = 1 / (cfs[j] * a))
  }
  results[[as.character(a)]] <- M
}

# Shared y-axis across all panels so the drop in CV with higher overdispersion
# is directly comparable.
ylim_all <- c(0, 1.02 * max(vapply(results, max, numeric(1))))
render <- function(dev_open) {
  dev_open()
  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
  for (a in alphas)
    boxplot(results[[as.character(a)]], names = cfs, ylim = ylim_all,
            main = paste("overdispersion", a),
            xlab = "correction factor", ylab = "CV", outline = FALSE)
  dev.off()
}
out <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "..", "output")
if (!dir.exists(out)) out <- "../output"
render(function() pdf(file.path(out, "pseudocount_simulation.pdf"), width = 9, height = 6))
render(function() png(file.path(out, "pseudocount_simulation.png"), width = 9, height = 6,
                      units = "in", res = 150))
cat("wrote", file.path(out, "pseudocount_simulation.{pdf,png}"), "\n")
