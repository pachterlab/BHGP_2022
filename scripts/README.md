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
