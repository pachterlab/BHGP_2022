# Figure 3: Angelidis pseudobulk PCA and PC1 loadings

This directory contains the source data and Python code for the Angelidis2019
pseudobulk analysis used in manuscript Figure 3.

Run:

```bash
python plot_angelidis_pseudobulk_figure.py
```

The script expects `scclr` and `edgePython` to be importable. It recomputes the
old-vs-young differential expression table from the raw pseudobulk counts.

The manuscript figure is written to:

```text
../../output/angelidis_pca_and_pc1_loadings_four_panel.pdf
```
