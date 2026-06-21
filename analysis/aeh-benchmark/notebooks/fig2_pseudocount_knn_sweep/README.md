# Mouse brain 10x input

The mouse-brain pseudocount/kNN sweep uses the public 10x Genomics
`neuron_1k_v3` filtered feature-barcode matrix:

```text
neuron_1k_v3_filtered_feature_bc_matrix.h5
```

The H5 matrix is not tracked in git. Download it with:

```bash
bash download_mouse_brain_neuron_1k_v3.sh
```

By default this writes to:

```text
../../benchmark/output/clr_local/data/mouse_brain_10x/neuron_1k_v3_filtered_feature_bc_matrix.h5
```

That location is under `benchmark/output/`, which is ignored by git. To use a
different location, pass it as the first argument:

```bash
bash download_mouse_brain_neuron_1k_v3.sh /path/to/data
```

Source URL:

```text
https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_filtered_feature_bc_matrix.h5
```

Expected SHA256:

```text
78166dda24a6103f6b690e0e5a392d89f9da429b84ee265ee4a0f4621f607695
```
