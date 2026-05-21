"""
Compare CLR-PCA vs standard logPF-PCA on scRNA-seq data (PBMC 3k).

Standard logPF-PCA: log(counts/total), center per gene, PCA
CLR-PCA:           log(counts/total), center per cell (Aitchison CLR), PCA
"""

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.sparse import issparse

sc.settings.verbosity = 0

# ── Load and filter ───────────────────────────────────────────────────────────
print("Loading PBMC 3k...")
adata = sc.datasets.pbmc3k()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Select top 2000 HVGs on log1p-normalized data (standard)
adata_norm = adata.copy()
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)
sc.pp.highly_variable_genes(adata_norm, n_top_genes=2000)
hvg_mask = adata_norm.var.highly_variable.values

# Raw counts for HVGs
counts = adata.X[:, hvg_mask]
if issparse(counts):
    counts = counts.toarray()
counts = counts.astype(float)
print(f"Shape after HVG selection: {counts.shape}")

# ── Proportional frequencies ──────────────────────────────────────────────────
totals = counts.sum(axis=1, keepdims=True)
PF = counts / totals

epsilon = 1e-8
logPF = np.log(PF + epsilon)   # cells x genes

# ── Two centering strategies ──────────────────────────────────────────────────
# Standard: center each gene across cells
logPF_gene_centered = logPF - logPF.mean(axis=0)

# CLR: center each cell across genes (Aitchison)
clr = logPF - logPF.mean(axis=1, keepdims=True)

# ── PCA ───────────────────────────────────────────────────────────────────────
n_components = 50
print("Running PCAs...")
pca_logpf = PCA(n_components=n_components, random_state=0).fit(logPF_gene_centered)
pca_clr   = PCA(n_components=n_components, random_state=0).fit(clr)

emb_logpf = pca_logpf.transform(logPF_gene_centered)
emb_clr   = pca_clr.transform(clr)

# ── Variance explained ────────────────────────────────────────────────────────
ev_logpf = pca_logpf.explained_variance_ratio_
ev_clr   = pca_clr.explained_variance_ratio_
print(f"\nVariance explained (top 10 PCs):")
print(f"  logPF-PCA: {100*ev_logpf[:10].sum():.1f}%  (PC1={100*ev_logpf[0]:.1f}%)")
print(f"  CLR-PCA:   {100*ev_clr[:10].sum():.1f}%  (PC1={100*ev_clr[0]:.1f}%)")

# ── PC1 correlation with sequencing depth ────────────────────────────────────
log_totals = np.log(totals.squeeze())
r_logpf = np.corrcoef(emb_logpf[:, 0], log_totals)[0, 1]
r_clr   = np.corrcoef(emb_clr[:, 0],   log_totals)[0, 1]
print(f"\nCorrelation of PC1 with log(total counts):")
print(f"  logPF-PCA: r = {r_logpf:.4f}")
print(f"  CLR-PCA:   r = {r_clr:.4f}")

# ── Clustering silhouette ─────────────────────────────────────────────────────
print("\nClustering (k=8):")
for emb, label in [(emb_logpf, "logPF-PCA"), (emb_clr, "CLR-PCA  ")]:
    km = KMeans(n_clusters=8, random_state=0, n_init=10)
    labels = km.fit_predict(emb[:, :20])
    sil = silhouette_score(emb[:, :20], labels, sample_size=2000, random_state=0)
    print(f"  {label}: silhouette = {sil:.4f}")

# ── Top genes in PC1 ─────────────────────────────────────────────────────────
gene_names = np.array(adata.var_names)[hvg_mask]
print(f"\nTop genes in PC1:")
for pca_model, label in [(pca_logpf, "logPF-PCA"), (pca_clr, "CLR-PCA  ")]:
    loadings = pca_model.components_[0]
    top = np.argsort(np.abs(loadings))[::-1][:5]
    print(f"  {label}: {[(gene_names[i], round(loadings[i],3)) for i in top]}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].plot(np.arange(1, 21), ev_logpf[:20] * 100, 'o-', label='logPF-PCA', color='steelblue')
axes[0].plot(np.arange(1, 21), ev_clr[:20]   * 100, 's-', label='CLR-PCA',   color='tomato')
axes[0].set_xlabel("PC"); axes[0].set_ylabel("Variance explained (%)")
axes[0].set_title("Scree plot"); axes[0].legend()

axes[1].scatter(log_totals, emb_logpf[:, 0], s=2, alpha=0.4,
                label=f"logPF-PCA  r={r_logpf:.3f}", color='steelblue')
axes[1].scatter(log_totals, emb_clr[:, 0],   s=2, alpha=0.4,
                label=f"CLR-PCA  r={r_clr:.3f}", color='tomato')
axes[1].set_xlabel("log(total counts)"); axes[1].set_ylabel("PC1 score")
axes[1].set_title("PC1 vs sequencing depth"); axes[1].legend(markerscale=4)

plt.tight_layout()
plt.savefig("clr_vs_logpf_pca.png", dpi=150, bbox_inches='tight')
print("\nSaved: clr_vs_logpf_pca.png")
