# PFlog / CLR benchmark and geometry analyses

This directory contains the adapted Ahlmann-Eltze and Huber transformation
benchmark with an added PFlog / shifted-CLR row, together with the additional
geometry analyses used in the manuscript and Supplementary Note.

- `benchmark/src/transformations/transformation_helper.R` defines `clr_fnc()`, the PFlog / CLR transform used by the benchmark framework.
- `benchmark/job_overview.yaml` adds `clr` to the consistency, simulation,
  downsampling, stratification, and siRNA recovery benchmark specifications.
- `benchmark/src/run_benchmarks.R` is the scheduler-oriented benchmark driver.
- `benchmark/run_clr_*.R` are local standalone runners used to fill CLR benchmark rows without running the full Slurm workflow.
- `notebooks/plot_benchmark_results.Rmd` renders the adapted Ahlmann-Eltze and
  Huber benchmark figure used as Supplementary Note Fig. 6.
- `notebooks/pseudocount_simulation.R` renders the pseudocount simulation used
  as Supplementary Note Fig. 2.
- `notebooks/fig2_pseudocount_knn_sweep/` contains the 10x Genomics mouse brain
  pseudocount and neighborhood-graph analysis used as Supplementary Note Fig. 3.
- `notebooks/fig3_angelidis_pseudobulk/` contains the main Fig. 2 Angelidis2019
  pseudobulk PCA, edgePython differential-expression, and PC1-loading comparison
  script and source data.
- `notebooks/fig4_lung_depth_geometry/` contains the main Fig. 3 lung
  cross-depth PC-loading analysis and the muscat 12-fold downsampling control
  quoted in the text.
- `notebooks/supplementary_figure_4/` contains the Seurat `"CLR"` depth-flaw
  plotting scripts and source summaries used as Supplementary Note Fig. 5.
- `scripts/plot_walk_depth_confounding.py` renders the linear-walk and
  random-walk depth-confounding diagnostic used as Supplementary Note Fig. 7.
- `extras/clr_vs_logpf_pca.py` is a standalone PBMC3k CLR-vs-logPF PCA
  diagnostic script retained for auditability.
