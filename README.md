# BHGP_2022

Code and analysis for *"Depth normalization for single-cell genomics count data"* (Booeshaghi, Hallgrímsdóttir, Gálvez-Merchán, Pachter).

## What's here

- `scripts/` — Python + bash pipeline: normalization (`norm_*.py`), per-method metrics (`metrics_*.py`), per-dataset plotting (`plot_*.py`), and k-NN consistency (`metrics_knn.py`).
- `analysis/` — Jupyter notebooks reproducing the per-dataset and aggregate figures.
- `data/datasets.txt` — list of GEO/SRA accessions used in the benchmark.
- `data/release_date.txt` — release date per accession.
- `pyproject.toml`, `uv.lock` — reproducible Python environment via [uv](https://docs.astral.sh/uv/).

## What's NOT here (intentionally)

The original matrix data (~1.3 TB of raw + normalized count matrices across 526 datasets) is not included to keep clones lightweight. Two ways to obtain it:

1. **Public sources**: each accession in `data/datasets.txt` is downloadable from GEO or SRA.
2. **Legacy data-bundled repo**: the original analysis with bundled intermediates is at [`pachterlab/BHGP_2022_v1`](https://github.com/pachterlab/BHGP_2022_v1).

## Quick start

```bash
git clone https://github.com/pachterlab/BHGP_2022
cd BHGP_2022
uv sync                                # install Python 3.11 + dependencies into .venv

# Place a raw matrix at data/<accession>/matrix.mtx.gz, then:
cd data && bash ../scripts/norm_sparse.sh
bash ../scripts/metrics_methods_sparse.sh
bash ../scripts/plot_sparse.sh
```

See `scripts/README.md` for the full 13-step pipeline order.

## License

BSD 2-Clause — see [`LICENSE`](./LICENSE).
