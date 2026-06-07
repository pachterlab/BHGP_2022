# Depth normalization for single-cell genomics count data

Code and data to reproduce *"Depth normalization for single-cell genomics count
data"* (Booeshaghi, Hallgrímsdóttir, Gálvez-Merchán, Pachter).

## Summary

Single-cell count normalization should **stabilize variance**, **remove
sequencing-depth effects**, and **preserve within-cell gene ranks**. We show that
a single transform satisfies all three (while being compatible with PCA and
invariant to gene order): proportional fitting, then `log1p`, then a second
proportional fitting — **PFlog1pPF**, which is exactly a **shifted centered
log-ratio (CLR)** transform. Across hundreds of datasets and the benchmarks of
[Ahlmann-Eltze & Huber (2023)](https://doi.org/10.1038/s41592-023-01814-1) it
matches or beats widely used alternatives (log1pPF, sctransform), most strikingly
in robustness to downsampling.

PFlog1pPF in a few lines (sparse in, dense out — the centering fills zeros):

```python
import numpy as np, scipy.sparse as sp
from scipy.io import mmread

def pf(mtx):                                 # proportional fitting to mean depth
    d = np.asarray(mtx.sum(1)).ravel()
    return sp.diags(d.mean() / d) @ mtx

mtx = mmread("raw.mtx.gz").tocsr()           # raw counts, cells x genes
log1ppf = pf(mtx); log1ppf.data = np.log1p(log1ppf.data)     # PF -> log1p
cell_mean = np.asarray(log1ppf.mean(1)).ravel()
pflog1ppf = log1ppf.toarray() - cell_mean[:, None]           # PF = per-cell centering
```

## Key results

- **PFlog1pPF (= shifted CLR)** is the unique count transform that is variance-
  stabilizing, depth-invariant, and rank-preserving (Supplementary Note).
- On **526 datasets** (437 passing QC) it is the best overall: it removes residual
  depth structure, keeps monotonicity, and stabilizes variance as well as log1pPF —
  whereas log1pPF/sctransform retain depth, and sctransform scrambles ranks.
- Reproducing the **Ahlmann-Eltze & Huber** k-NN benchmarks, PFlog1pPF is strongest
  on downsampling: ~**36.8/50** neighbors recovered vs ~5.8 for other methods.
- Its embedding geometry is **more stable to feature-panel choice** than log1pPF.
- Seurat's `"CLR"` is **not** the centered log-ratio and does not remove depth.

## Repository layout

- **`scripts/`** — Python/bash pipeline for the 526-dataset analysis (Fig. 1):
  normalization (`norm_*`), per-method and cell-type metrics (`metrics_*`), and
  plotting (`plot_*`). Pipeline order documented in [`scripts/README.md`](scripts/README.md).
- **`PFlogPF/`** — R reproduction of the Ahlmann-Eltze & Huber k-NN benchmark
  (Fig. 2) and the downstream/supplementary figures (downsampling, feature-panel
  stability, Seurat-`"CLR"`, pseudocount simulation). `PFlogPF/benchmark/` holds the
  benchmark framework with committed result tables and a pinned `renv.lock`.
- **`analysis/`** — notebooks and the rendered manuscript figures in `analysis/figures/`.
- **`data/`** — `datasets.txt` (526 GEO/SRA accessions), `release_date.txt`, and
  matrix-size tables.

## Data

The raw and normalized count matrices (~1 TB across 526 datasets) are not bundled.
Each accession in `data/datasets.txt` is available from GEO/SRA; the original
analysis with bundled intermediates is archived at
[`pachterlab/BHGP_2022_v1`](https://github.com/pachterlab/BHGP_2022_v1).

## Reproducing

```bash
uv sync                 # Python 3.11 env into .venv (see pyproject.toml / uv.lock)
```

- **Fig. 1 (526-dataset benchmark):** run the `scripts/` pipeline per
  `scripts/README.md` (normalize → metrics → `plot_summary.py` / `plot_bar.py`).
- **Fig. 2 (k-NN benchmark) + supplement:** R, under `PFlogPF/` (restore `renv`,
  run `benchmark/`, render `notebooks/plot_benchmark_results.Rmd`).

## License

BSD 2-Clause — see [`LICENSE`](./LICENSE).
