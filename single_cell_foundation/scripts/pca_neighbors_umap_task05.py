import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

print("Task 05: PCA, nearest-neighbor graph, and UMAP")

# -----------------------------
# 1. Paths
# -----------------------------
input_h5ad = Path("data/pbmc3k_normalized_hvg.h5ad")
data_dir = Path("data")
results_dir = Path("results/task05")
figures_dir = Path("figures/task05")

results_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

if not input_h5ad.exists():
    raise FileNotFoundError(f"Missing input file: {input_h5ad}")

# -----------------------------
# 2. Load AnnData
# -----------------------------
adata = sc.read_h5ad(input_h5ad)

print("\nInput AnnData:")
print(adata)

n_cells, n_genes = adata.shape
n_hvg = int(adata.var["highly_variable"].sum())

print(f"Cells: {n_cells}")
print(f"Genes: {n_genes}")
print(f"Highly variable genes: {n_hvg}")

# -----------------------------
# 3. Scale expression values
# -----------------------------
# PCA works better when genes are put on comparable scales.
# This changes adata.X, so log-normalized values remain preserved in adata.raw.
adata.raw = adata

sc.pp.scale(adata, max_value=10)

# -----------------------------
# 4. PCA using highly variable genes
# -----------------------------
sc.tl.pca(
    adata,
    n_comps=50,
    svd_solver="arpack",
    mask_var="highly_variable"
)

# -----------------------------
# 5. Nearest-neighbor graph
# -----------------------------
n_neighbors = 10
n_pcs = 40

sc.pp.neighbors(
    adata,
    n_neighbors=n_neighbors,
    n_pcs=n_pcs
)

# -----------------------------
# 6. UMAP
# -----------------------------
sc.tl.umap(adata)

# -----------------------------
# 7. Save processed AnnData
# -----------------------------
output_h5ad = data_dir / "pbmc3k_pca_neighbors_umap.h5ad"
adata.write(output_h5ad)
print(f"\nSaved AnnData to: {output_h5ad}")

# -----------------------------
# 8. Save summary
# -----------------------------
summary = pd.DataFrame({
    "item": [
        "cells",
        "genes",
        "highly_variable_genes",
        "pca_n_comps",
        "neighbors_n_neighbors",
        "neighbors_n_pcs",
        "umap_dimensions"
    ],
    "value": [
        n_cells,
        n_genes,
        n_hvg,
        50,
        n_neighbors,
        n_pcs,
        2
    ]
})

summary.to_csv(results_dir / "task05_pca_neighbors_umap_summary.tsv", sep="\t", index=False)

# PCA variance ratio
pca_variance = pd.DataFrame({
    "PC": [f"PC{i+1}" for i in range(len(adata.uns["pca"]["variance_ratio"]))],
    "variance_ratio": adata.uns["pca"]["variance_ratio"],
    "variance": adata.uns["pca"]["variance"]
})

pca_variance.to_csv(results_dir / "pca_variance_ratio.tsv", sep="\t", index=False)

# UMAP coordinates
umap_df = pd.DataFrame(
    adata.obsm["X_umap"],
    index=adata.obs_names,
    columns=["UMAP1", "UMAP2"]
)

umap_df.to_csv(results_dir / "umap_coordinates.tsv", sep="\t")

# PCA coordinates, first 10 PCs
pca_df = pd.DataFrame(
    adata.obsm["X_pca"][:, :10],
    index=adata.obs_names,
    columns=[f"PC{i+1}" for i in range(10)]
)

pca_df.to_csv(results_dir / "pca_coordinates_first10.tsv", sep="\t")

# -----------------------------
# 9. Figures
# -----------------------------

# PCA variance ratio plot
plt.figure(figsize=(7, 4))
plt.plot(
    range(1, 51),
    adata.uns["pca"]["variance_ratio"][:50],
    marker="o",
    markersize=3
)
plt.xlabel("Principal component")
plt.ylabel("Variance ratio")
plt.title("PCA variance explained")
plt.tight_layout()
plt.savefig(figures_dir / "pca_variance_ratio.png", dpi=150)
plt.close()

# UMAP colored by QC metrics
for color in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
    sc.pl.umap(
        adata,
        color=color,
        show=False
    )
    plt.savefig(figures_dir / f"umap_{color}.png", dpi=150, bbox_inches="tight")
    plt.close()

# UMAP without color
plt.figure(figsize=(5, 4))
plt.scatter(
    adata.obsm["X_umap"][:, 0],
    adata.obsm["X_umap"][:, 1],
    s=5,
    alpha=0.7
)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("PBMC3k UMAP")
plt.tight_layout()
plt.savefig(figures_dir / "umap_plain.png", dpi=150)
plt.close()

print("\nSummary:")
print(summary)

print("\nMain outputs:")
print(f"  {output_h5ad}")
print(f"  {results_dir / 'task05_pca_neighbors_umap_summary.tsv'}")
print(f"  {results_dir / 'pca_variance_ratio.tsv'}")
print(f"  {results_dir / 'umap_coordinates.tsv'}")
print(f"  {figures_dir}")
print("Task 05 finished.")
