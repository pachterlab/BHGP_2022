# PFlogPF / CLR Benchmark Additions

This directory contains the code and manuscript figure outputs for the Ahlmann-Eltze and Huber transformation benchmark plots with the PFlogPF / CLR transformation added.

- `benchmark/src/transformations/transformation_helper.R` defines `clr_fnc()`, the PFlogPF / CLR transform used by the benchmark framework.
- `benchmark/job_overview.yaml` adds `clr` to the consistency, simulation, downsampling, stratification, and siRNA recovery benchmark specifications.
- `benchmark/src/run_benchmarks.R` is the scheduler-oriented benchmark driver.
- `benchmark/run_clr_*.R` are local standalone runners used to fill CLR benchmark rows without running the full Slurm workflow.
- `benchmark/src/downsampling_benchmark/calculate_sirna_quantitative_recovery.R` and related files implement the siRNA quantitative recovery analysis.
- `notebooks/annotation_helper.R` and the Rmd/R plotting scripts add CLR labels and render the benchmark figures.
- `notebooks/fig2_with_clr/` contains the compact figure-specific reproduction scripts.
- `notebooks/fig3_panel_sensitivity/` contains the Figure 3 Smart-seq3 feature-panel sensitivity source data and scatter-panel plotting script.
- `notebooks/supplementary_figure_4/` contains the Supplementary Figure 4 Seurat CLR depth-flaw plotting script and source summaries.
- `output/smartseq3_fibroblasts_paired_procrustes_displacement_projection.pdf` and `output/smartseq3_fibroblasts_paired_distance_scatter.pdf` are the Figure 3 panels used in the manuscript.
- `extras/clr_vs_logpf_pca.py` is a standalone PBMC3k CLR-vs-logPF PCA sanity-check script.
