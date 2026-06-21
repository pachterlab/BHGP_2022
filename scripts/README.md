# 526-dataset normalization pipeline

This directory contains the Python, R, and shell scripts used for the
526-dataset normalization benchmark summarized in main Fig. 1. The scripts assume
one directory per observation or dataset, denoted below as `OBS`.

## Per-dataset processing

Prepare the sparse raw-count input:

```bash
mkdir -p OBS/sparse
cp OBS/matrix.mtx.gz OBS/sparse/raw.mtx.gz
```

Then run the sparse-matrix metrics and normalizations:

1. From `OBS/sparse`, `../../scripts/metrics_matrix.sh` computes raw matrix
   metrics.
2. From `OBS/sparse`, `../../scripts/norm_sparse.sh` computes sparse
   normalizations.
3. From `OBS/sparse`, `../../scripts/metrics_methods_sparse.sh` computes
   per-method metrics for sparse normalizations.

Run sctransform and the dense-method metrics:

4. From `OBS`, `../scripts/norm_sctransform.sh` computes the sctransform matrix
   and writes the raw matrix used for the sctransform-selected gene subset.
5. From `OBS`, `../scripts/norm_sparse.sh` computes sparse normalizations on
   that gene subset.
6. From `OBS`, `../scripts/norm_cp10k_log_scale.sh` computes the log1p(CP10K)
   matrix with scaling.
7. From `OBS`, `../scripts/metrics_matrix.sh` computes matrix metrics on the
   subset.
8. From `OBS`, `../scripts/metrics_methods_sparse.sh` computes sparse-method
   metrics on the subset.
9. From `OBS`, `../scripts/metrics_methods_dense.sh` computes dense-method
   metrics for sctransform and log1p(CP10K) with scaling.

Plot per-dataset diagnostics:

10. From `OBS`, `../scripts/plot_sparse.sh` renders sparse-method comparison
    plots.
11. From `OBS`, `../scripts/plot_all.sh` renders comparison plots over all
    methods.

## Manuscript figures

The per-dataset processing produces
`OBS/subset_genes/OBS_subset_genes_metrics.json`, which is the input for the
main Fig. 1 summary plots. Let `DATA` denote the root containing all dataset
directories.

Main Fig. 1a, the Angelidis2019 Type 2 pneumocyte cell-type metrics:

```bash
SC_CLR_PATH=/path/to/scclr/python \
python scripts/metrics_celltype.py DATA/angelidis_2019 Type_2_pneumocytes \
  analysis/figures/angelidis_celltype_metrics.json \
  --clr_calibrate --scclr-path /path/to/scclr/python

python scripts/plot_celltype_bar.py \
  analysis/figures/angelidis_celltype_metrics.json \
  analysis/figures/angelidis_celltype_bar
```

Main Fig. 1b, the across-dataset summary:

```bash
python scripts/plot_summary.py DATA analysis/figures/summary_subset_genes angelidis_2019
```

PFlog rows are computed with the current `scclr` shifted CLR transform, stored
sparsely as `center(log1p(4 * alpha * x))`. This is algebraically equivalent to
centering `log(x + 1/(4 * alpha))` because the omitted `log(4 * alpha)` constant
cancels under within-cell centering.

Legacy `pf_log_pf.mtx.gz` files are retained only as historical artifacts and
are not used for the manuscript PFlog results.

`scripts/plot_bar.py` renders an all-cells single-dataset metric bar
`OBS_metrics-bar` from a `*_subset_genes_metrics.json`; this is a per-dataset
diagnostic, not a numbered main figure.
