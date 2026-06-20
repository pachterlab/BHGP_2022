OBS is the "id" for an observation

NB: for points 3, 4, and 5 the path in the file will need to be `$1/sparse/` and will need to be changed to `$1/` for points 7, 9, and 10.

1. copy over matrix.mtx.gz and kneeplot.png to OBS
2. make sparse folder for each observation `mkdir -p OBS/sparse`, `cp OBS/matrix.mtx.gz OBS/sparse/raw.mtx.gz`
3. `../scripts/metrics_matrix.sh` - computes matrix metrics for sparse dudes
4. `../scripts/norm_sparse.sh` - computes all normalizations on sparse dudes
5. `../scripts/metrics_methods_sparse.sh` - computes method metrics for sparse dudes
6. `../scripts/norm_sctransform.sh` - compute sctransform (also produces raw) (gzip raw)
7. `../scripts/norm_sparse.sh` - compute sparse normalization for raw (that was produced in the previous step)
8. `../scripts/norm_cp10k_log_scale.sh` - compute cp10k log scale
9. `../scripts/metrics_matrix.sh` - compute matrix metrics on raw (from sct) 
10. `../scripts/metrics_methods_sparse.sh` - compute method metrics sparse dudes
11. `../scripts/metrics_methods_dense.sh` - compute method metrics for sctransform and cp10klogscale
12. `../scripts/plot_sparse.sh` - produce method comparison plot (on sparse only)
13. `..scripts/plot_all.sh` - produce method comparison plot (on all)

## Manuscript figures

Steps 1-13 produce, per dataset, `<DS>/subset_genes/<DS>_subset_genes_metrics.json`
(the per-method metric inputs). The manuscript figures are then rendered from
those JSONs with the scripts below (DATA = the data root containing the dataset
dirs, e.g. `synchromesh/data`):

- **Fig 1a** — angelidis_2019 within-cell-type per-method bars (PC1-loading
  entropy / # false-positive DE genes / 1-|within-celltype Spearman|) for the
  Type_2_pneumocytes cell type. Compute the metrics JSON, then plot:
  `SC_CLR_PATH=/path/to/scclr/python python scripts/metrics_celltype.py DATA/angelidis_2019 Type_2_pneumocytes analysis/figures/angelidis_celltype_metrics.json --clr_calibrate --scclr-path /path/to/scclr/python`
  `python scripts/plot_celltype_bar.py analysis/figures/angelidis_celltype_metrics.json analysis/figures/angelidis_celltype_bar`
- **Fig 1b** — across-dataset summary boxplot over all datasets' metrics JSONs:
  `python scripts/plot_summary.py DATA analysis/figures/summary_subset_genes angelidis_2019`

PFlog rows are computed with the corrected runorm/scclr transform,
`center(log1p(4*alpha*x))`. Legacy `pf_log_pf.mtx.gz` files are kept only as
historical artifacts and should not be used as the manuscript PFlog result.

(`scripts/plot_bar.py` renders an all-cells single-dataset metric bar
`<DS>_metrics-bar` from a `*_subset_genes_metrics.json` — a per-dataset analog of
Fig 1b, not a numbered main figure.)
