import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

print("Task 06: Leiden clustering")

# -----------------------------
# 1. Paths
# -----------------------------
input_h5ad = Path("data/pbmc3k_pca_neighbors_umap.h5ad")
data_dir = Path("data")
results_dir = Path("results/task06")
figures_dir = Path("figures/task06")

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

if "X_umap" not in adata.obsm:
    raise ValueError("UMAP coordinates not found. Run Task 05 first.")

if "neighbors" not in adata.uns:
    raise ValueError("Neighbor graph not found. Run Task 05 first.")

# -----------------------------
# 3. Leiden clustering at multiple resolutions
# -----------------------------
# Resolution controls clustering granularity.
# Lower resolution: fewer, larger clusters.
# Higher resolution: more, smaller clusters.
resolutions = [0.3, 0.5, 0.8, 1.0]

for res in resolutions:
    key = f"leiden_r{str(res).replace('.', '_')}"
    sc.tl.leiden(
        adata,
        resolution=res,
        key_added=key,
        random_state=0,
        flavor="igraph",
        n_iterations=2,
        directed=False
    )
    print(f"{key}: {adata.obs[key].nunique()} clusters")

# Use resolution 0.8 as the primary beginner-level clustering.
adata.obs["leiden"] = adata.obs["leiden_r0_8"].copy()

n_clusters_primary = adata.obs["leiden"].nunique()
print(f"\nPrimary Leiden clusters: {n_clusters_primary}")

# -----------------------------
# 4. Save cluster count tables
# -----------------------------
cluster_summary_rows = []

for res in resolutions:
    key = f"leiden_r{str(res).replace('.', '_')}"
    counts = adata.obs[key].value_counts().sort_index()
    counts_df = counts.rename_axis("cluster").reset_index(name="n_cells")
    counts_df["resolution"] = res
    counts_df.to_csv(results_dir / f"{key}_cluster_counts.tsv", sep="\t", index=False)

    cluster_summary_rows.append({
        "resolution": res,
        "key": key,
        "n_clusters": adata.obs[key].nunique(),
        "smallest_cluster_size": int(counts.min()),
        "largest_cluster_size": int(counts.max())
    })

cluster_summary = pd.DataFrame(cluster_summary_rows)
cluster_summary.to_csv(results_dir / "leiden_resolution_summary.tsv", sep="\t", index=False)

primary_counts = adata.obs["leiden"].value_counts().sort_index()
primary_counts_df = primary_counts.rename_axis("cluster").reset_index(name="n_cells")
primary_counts_df.to_csv(results_dir / "leiden_primary_cluster_counts.tsv", sep="\t", index=False)

# Save cell-to-cluster table
cell_cluster_df = adata.obs[
    ["leiden", "leiden_r0_3", "leiden_r0_5", "leiden_r0_8", "leiden_r1_0"]
].copy()

cell_cluster_df.to_csv(results_dir / "cell_leiden_clusters.tsv", sep="\t")

# -----------------------------
# 5. Save summary
# -----------------------------
summary = pd.DataFrame({
    "item": [
        "cells",
        "primary_resolution",
        "primary_cluster_key",
        "primary_n_clusters",
        "tested_resolutions"
    ],
    "value": [
        adata.n_obs,
        0.8,
        "leiden",
        n_clusters_primary,
        ",".join(map(str, resolutions))
    ]
})

summary.to_csv(results_dir / "task06_leiden_clustering_summary.tsv", sep="\t", index=False)

# -----------------------------
# 6. Figures
# -----------------------------

# UMAP colored by primary Leiden clusters
sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    title="Leiden clusters, resolution 0.8",
    show=False
)
plt.savefig(figures_dir / "umap_leiden_primary_r0_8.png", dpi=150, bbox_inches="tight")
plt.close()

# UMAP colored by different resolutions
for res in resolutions:
    key = f"leiden_r{str(res).replace('.', '_')}"
    sc.pl.umap(
        adata,
        color=key,
        legend_loc="on data",
        title=f"Leiden clusters, resolution {res}",
        show=False
    )
    plt.savefig(figures_dir / f"umap_{key}.png", dpi=150, bbox_inches="tight")
    plt.close()

# Cluster size barplot for primary clustering
plt.figure(figsize=(7, 4))
plt.bar(primary_counts_df["cluster"].astype(str), primary_counts_df["n_cells"])
plt.xlabel("Leiden cluster")
plt.ylabel("Number of cells")
plt.title("Primary Leiden cluster sizes, resolution 0.8")
plt.tight_layout()
plt.savefig(figures_dir / "leiden_primary_cluster_sizes.png", dpi=150)
plt.close()

# -----------------------------
# 7. Save AnnData
# -----------------------------
output_h5ad = data_dir / "pbmc3k_leiden.h5ad"
adata.write(output_h5ad)

print("\nSummary:")
print(summary)

print("\nPrimary cluster counts:")
print(primary_counts_df)

print("\nResolution summary:")
print(cluster_summary)

print("\nMain outputs:")
print(f"  {output_h5ad}")
print(f"  {results_dir / 'task06_leiden_clustering_summary.tsv'}")
print(f"  {results_dir / 'leiden_primary_cluster_counts.tsv'}")
print(f"  {results_dir / 'leiden_resolution_summary.tsv'}")
print(f"  {figures_dir}")
print("Task 06 finished.")
