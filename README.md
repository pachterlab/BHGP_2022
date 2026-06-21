# Normalization for negative-binomial count data

Code and data to reproduce *"Normalization for negative-binomial count
data"* (Booeshaghi, Hallgrímsdóttir, Gálvez-Merchán, Pachter).

## Summary

Single-cell count normalization should **stabilize technical variance**,
**remove sequencing-depth effects**, and **preserve within-cell gene ranks**.
PFlog addresses these requirements by applying an Anscombe-calibrated shifted
logarithm followed by within-cell CLR centering. Equivalently, it is a **shifted
centered log-ratio (CLR)** transform. Across hundreds of datasets, and in an
adapted analysis of the benchmarks of
[Ahlmann-Eltze & Huber (2023)](https://doi.org/10.1038/s41592-023-01814-1),
PFlog removes residual depth structure while preserving rank information and
stabilizing technical variance.

The shift is set by the **delta method**. For a negative-binomial mean-variance
model `var = μ + α·μ²`, the Anscombe-derived count-scale pseudocount is
`y₀ = 1/(4α)`. In software notation, for cell `c` with depth `s_c`, this is
equivalent to the cell-specific normalized shifted-log scale
`K_c = 4·α·s_c`. At a representative depth `s_*`, the corresponding plotted
scale is `K_* = 4·α·s_*`.

PFlog in a few lines, written in its count-scale pseudocount form:

```python
import numpy as np
from scipy.io import mmread
from scipy.optimize import curve_fit

def overdispersion(mtx):                     # fit var = mu + alpha*mu^2 across genes
    mu  = np.asarray(mtx.mean(0)).ravel()
    m2  = mtx.copy(); m2.data **= 2
    var = np.asarray(m2.mean(0)).ravel() - mu**2
    return curve_fit(lambda x, a: x + a * x**2, mu, var)[0][0]

mtx = mmread("raw.mtx.gz").tocsr()           # raw counts, cells x genes
alpha = overdispersion(mtx)
y0 = 1 / (4 * alpha)

shifted_log = np.log(mtx.toarray() + y0)
cell_mean = shifted_log.mean(axis=1)
pflog = shifted_log - cell_mean[:, None]
```

The sparse implementation used in the manuscript computes the equivalent
`center(log1p(4 * alpha * x))` form and stores the shifted-log sparse matrix
together with its centering vector, avoiding materialization of the dense
centered matrix when possible. For a cell-type-specific analysis, estimate
`α` within the cell type; pooled-dataset estimates can be inflated by
between-cell-type heterogeneity.

## Key results

- **PFlog (= shifted CLR)** combines an Anscombe-calibrated count-scale
  pseudocount with CLR centering, giving technical variance stabilization,
  depth invariance, and rank preservation.
- On **526 datasets** (437 passing QC) it is the best overall: it removes residual
  depth structure, keeps monotonicity, and stabilizes variance as well as log1pPF —
  whereas log1pPF/sctransform retain depth, and sctransform scrambles ranks.
- In Angelidis2019 pseudobulk profiles, both the count-scale pseudocount and CLR
  centering materially improve the agreement between unsupervised PC1 loadings
  and supervised old-vs-young differential expression.
- In lung datasets with large depth differences, PFlog improves cross-assay PC1
  loading agreement relative to fixed-scale log1p(CP10K) and to shifted log
  without CLR centering.
- In an adapted **Ahlmann-Eltze & Huber** benchmark, PFlog is added as a
  held-out CLR/PFlog row and performs favorably on retained simulation and
  downsampling comparisons.
- Seurat's `"CLR"` is **not** the centered log-ratio and does not remove depth.

## Repository layout

- **`scripts/`** — Python/bash pipeline for the 526-dataset analysis (Fig. 1):
  normalization (`norm_*`), per-method and cell-type metrics (`metrics_*`), and
  plotting (`plot_*`). Pipeline order documented in [`scripts/README.md`](scripts/README.md).
- **`analysis/aeh-benchmark/`** — adapted reproduction of the Ahlmann-Eltze & Huber
  k-NN benchmark and the manuscript/supplement analyses for Angelidis2019
  pseudobulk PCA, lung cross-depth loading stability, mouse-brain pseudocount
  sweeps, Seurat-`"CLR"`, pseudocount simulation, and walk-simulation depth
  diagnostics. `analysis/aeh-benchmark/benchmark/` holds the
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
- **Main Fig. 2, main Fig. 3, and Supplementary Note figures:** use the
  figure-specific scripts under `analysis/aeh-benchmark/notebooks/` and the
  adapted benchmark workflow under `analysis/aeh-benchmark/benchmark/`.

## License

BSD 2-Clause — see [`LICENSE`](./LICENSE).
