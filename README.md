# Depth normalization for single-cell genomics count data

Code and data to reproduce *"Depth normalization for single-cell genomics count
data"* (Booeshaghi, Hallgrímsdóttir, Gálvez-Merchán, Pachter).

## Summary

Single-cell count normalization should **stabilize variance**, **remove
sequencing-depth effects**, and **preserve within-cell gene ranks**. We show that
a single transform satisfies all three (while being compatible with PCA and
invariant to gene order): proportional fitting, then `log1p`, then centering —
**PFlog**, which is exactly a **shifted centered log-ratio (CLR)** transform. Across hundreds of datasets and the benchmarks of
[Ahlmann-Eltze & Huber (2023)](https://doi.org/10.1038/s41592-023-01814-1) it
matches or beats widely used alternatives (log1pPF, sctransform), most strikingly
in robustness to downsampling.

The shift (pseudocount) is set by the **delta method**: the first proportional
fitting targets `K = 4·α·s`, where `α` is the dataset's negative-binomial
overdispersion (`var = μ + α·μ²`) and `s` is the mean total UMI per cell. This is
exactly the variance-stabilizing count-scale pseudocount `y₀ = 1/(4α)` (see the
Supplementary Note); larger overdispersion ⇒ smaller pseudocount.

PFlog in a few lines (sparse in, dense out — the centering fills zeros):

```python
import numpy as np, scipy.sparse as sp
from scipy.io import mmread
from scipy.optimize import curve_fit

def overdispersion(mtx):                     # fit var = mu + alpha*mu^2 across genes
    mu  = np.asarray(mtx.mean(0)).ravel()
    m2  = mtx.copy(); m2.data **= 2
    var = np.asarray(m2.mean(0)).ravel() - mu**2
    return curve_fit(lambda x, a: x + a * x**2, mu, var)[0][0]

def pf(mtx, target):                         # proportional fitting to `target` total
    d = np.asarray(mtx.sum(1)).ravel()
    return sp.diags(target / d) @ mtx

mtx = mmread("raw.mtx.gz").tocsr()           # raw counts, cells x genes
s   = mtx.sum(1).mean()                      # scale factor: mean total UMI / cell
K   = 4 * overdispersion(mtx) * s            # delta-method first-PF constant (y0 = 1/(4a))

log1ppf = pf(mtx, K); log1ppf.data = np.log1p(log1ppf.data)  # PF to K -> log1p
cell_mean = np.asarray(log1ppf.mean(1)).ravel()
pflog = log1ppf.toarray() - cell_mean[:, None]               # per-cell centering
```

For a cell-type-specific analysis, estimate `α` and `s` *within* the cell type —
the pooled-dataset `α` is inflated by between-cell-type heterogeneity.

## Key results

- **PFlog (= shifted CLR)** is the unique count transform that is variance-
  stabilizing, depth-invariant, and rank-preserving (Supplementary Note).
- On **526 datasets** (437 passing QC) it is the best overall: it removes residual
  depth structure, keeps monotonicity, and stabilizes variance as well as log1pPF —
  whereas log1pPF/sctransform retain depth, and sctransform scrambles ranks.
- Reproducing the **Ahlmann-Eltze & Huber** k-NN benchmarks, PFlog is strongest
  on downsampling: ~**36.8/50** neighbors recovered vs ~5.8 for other methods.
- Its embedding geometry is **more stable to feature-panel choice** than log1pPF.
- Seurat's `"CLR"` is **not** the centered log-ratio and does not remove depth.

## Repository layout

- **`scripts/`** — Python/bash pipeline for the 526-dataset analysis (Fig. 1):
  normalization (`norm_*`), per-method and cell-type metrics (`metrics_*`), and
  plotting (`plot_*`). Pipeline order documented in [`scripts/README.md`](scripts/README.md).
- **`analysis/aeh-benchmark/`** — R reproduction of the Ahlmann-Eltze & Huber k-NN benchmark
  (Fig. 2) and the downstream/supplementary figures (downsampling, feature-panel
  stability, Seurat-`"CLR"`, pseudocount simulation). `analysis/aeh-benchmark/benchmark/` holds the
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
- **Fig. 2 (k-NN benchmark) + supplement:** R, under `analysis/aeh-benchmark/` (restore `renv`,
  run `benchmark/`, render `notebooks/plot_benchmark_results.Rmd`).

## License

BSD 2-Clause — see [`LICENSE`](./LICENSE).
