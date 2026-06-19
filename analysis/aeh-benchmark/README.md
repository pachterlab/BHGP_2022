# PFlog / CLR Benchmark Additions

This directory contains the code and manuscript figure outputs for the Ahlmann-Eltze and Huber transformation benchmark plots with the PFlog / CLR transformation added.

- `benchmark/src/transformations/transformation_helper.R` defines `clr_fnc()`, the PFlog / CLR transform used by the benchmark framework.
- `benchmark/job_overview.yaml` adds `clr` to the consistency, simulation, downsampling, stratification, and siRNA recovery benchmark specifications.
- `benchmark/src/run_benchmarks.R` is the scheduler-oriented benchmark driver.
- `benchmark/run_clr_*.R` are local standalone runners used to fill CLR benchmark rows without running the full Slurm workflow.
- `benchmark/src/downsampling_benchmark/calculate_sirna_quantitative_recovery.R` and related files implement the siRNA quantitative recovery analysis.
- `notebooks/annotation_helper.R` and the Rmd/R plotting scripts add CLR labels and render the benchmark figures.
- `notebooks/fig2_with_clr/` contains the compact figure-specific reproduction scripts.
- `notebooks/fig3_angelidis_pseudobulk/` contains the manuscript Figure 3 Angelidis2019 pseudobulk PCA, edgePython differential-expression, and PC1-loading comparison script and source data.
- `notebooks/supplementary_figure_4/` contains the Supplementary Note Figure 6 Seurat CLR depth-flaw plotting scripts and source summaries.
- `notebooks/supplementary_figure_6/` contains the Supplementary Note Figure 7 siRNA PCA-dimensionality plotting script and source summary.
- `output/angelidis_pca_and_pc1_loadings_four_panel.pdf` is Figure 3 in the manuscript.
- `extras/clr_vs_logpf_pca.py` is a standalone PBMC3k CLR-vs-logPF PCA sanity-check script.
