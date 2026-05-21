# PFlogPF / CLR Benchmark Additions

This directory contains the code needed to reproduce the Ahlmann-Eltze and Huber transformation benchmark plots with the PFlogPF / CLR transformation added.

The code was copied from the local `transformGamPoi-Paper` analysis tree and keeps the same relative layout where possible:

- `benchmark/src/transformations/transformation_helper.R` defines `clr_fnc()`, the PFlogPF / CLR transform used by the benchmark framework.
- `benchmark/job_overview.yaml` adds `clr` to the consistency, simulation, downsampling, stratification, and siRNA recovery benchmark specifications.
- `benchmark/src/run_benchmarks.R` is the scheduler-oriented benchmark driver.
- `benchmark/run_clr_*.R` are local standalone runners used to fill CLR benchmark rows without running the full Slurm workflow.
- `benchmark/src/downsampling_benchmark/calculate_sirna_quantitative_recovery.R` and related files implement the siRNA quantitative recovery analysis.
- `notebooks/annotation_helper.R` and the copied Rmd/R plotting scripts add CLR labels and render the benchmark figures.
- `notebooks/fig2_with_clr/` contains the compact figure-specific reproduction scripts.
- `extras/clr_vs_logpf_pca.py` is a standalone PBMC3k CLR-vs-logPF PCA sanity-check script.

Large generated artifacts are intentionally not included here: downloaded datasets, benchmark logs, intermediate result stores, rendered PDFs/PNGs/HTML, and vendored `renv/library` package files.

