# Main Fig. 3: lung cross-depth PC loading geometry

This directory contains the source data and Python code for main Fig. 3 and the
muscat 12-fold downsampling kNN result quoted in the Discussion.

## Inputs

- `data/hlca_lung_seqwell_10xv3_balanced_raw.h5ad`: raw-count HLCA lung subset
  with 560 Seq-Well cells and 560 10x 3' v3 cells, balanced across eight cell
  types.
- `data/muscat_seed1_counts.mtx`: raw muscat seed-1 count matrix used for the
  12x downsampling control.
- `data/muscat_seed1_truth_cluster_labels.tsv`: coarse labels derived from the
  muscat truth matrix for the supervised scANVI step.

## Main Fig. 3

By default, the HLCA PC loading comparison uses the top 1,000 common genes by
combined mean expression after requiring detection in more than 1% of cells in
both assays.

Run from the repository root:

```bash
python analysis/aeh-benchmark/notebooks/fig4_lung_depth_geometry/plot_hlca_lung_pc1_loading_stability.py
```

If `scclr` is not installed, pass a local checkout:

```bash
python analysis/aeh-benchmark/notebooks/fig4_lung_depth_geometry/plot_hlca_lung_pc1_loading_stability.py \
  --scclr-path /path/to/scclr/python
```

Outputs:

- `figures/lung_10x_vs_seqwell_depth_stability_pc1_loadings.pdf`
- `figures/lung_10x_vs_seqwell_depth_stability_pc1_loadings.png`
- `tables/lung_10x_vs_seqwell_depth_stability_summary.tsv`
- `tables/lung_10x_vs_seqwell_depth_stability_pc1_loadings.tsv`

The manuscript copy of the figure is also stored at:

- `../../output/lung_10x_vs_seqwell_depth_stability_pc1_loadings.pdf`

## muscat k=10 downsampling control

Run from the repository root in an environment with `scvi-tools`, `torch`,
`anndata`, `sklearn`, and `scclr`:

```bash
python analysis/aeh-benchmark/notebooks/fig4_lung_depth_geometry/run_muscat_downsample_k10_recovery.py
```

If `scclr` is not installed, pass a local checkout as above. The output table is:

- `tables/muscat_downsample12_k10_recovery.tsv`

The key manuscript values are the `full_knn_recovery_from_downsampled` entries
for `PFlog PCA` (`0.3100`) and `scANVI latent` (`0.1716`).
