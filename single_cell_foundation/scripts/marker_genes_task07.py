import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

print("Task 07: marker gene detection")

# -----------------------------
# 1. Paths
# -----------------------------
input_h5ad = Path("data/pbmc3k_leiden.h5ad")
data_dir = Path("data")
results_dir = Path("results/task07")
figures_dir = Path("figures/task07")

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

if "leiden" not in adata.obs.columns:
    raise ValueError("Missing adata.obs['leiden']. Run Task 06 first.")

if adata.raw is None:
    raise ValueError("Missing adata.raw. Task 05 should preserve log-normalized values in adata.raw.")

print("\nCluster counts:")
print(adata.obs["leiden"].value_counts().sort_index())

# -----------------------------
# 3. Detect marker genes
# -----------------------------
# Use adata.raw because adata.X was scaled in Task 05.
# adata.raw contains log-normalized expression values before scaling.
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    use_raw=True,
    key_added="rank_genes_leiden"
)

# -----------------------------
# 4. Convert marker results to tables
# -----------------------------
rg = adata.uns["rank_genes_leiden"]
groups = rg["names"].dtype.names

all_marker_tables = []

for group in groups:
    df = pd.DataFrame({
        "cluster": group,
        "gene": rg["names"][group],
        "score": rg["scores"][group],
        "logfoldchange": rg["logfoldchanges"][group],
        "pval": rg["pvals"][group],
        "pval_adj": rg["pvals_adj"][group],
    })

    df["rank"] = range(1, len(df) + 1)

    # Reorder columns
    df = df[
        ["cluster", "rank", "gene", "score", "logfoldchange", "pval", "pval_adj"]
    ]

    df.to_csv(results_dir / f"cluster_{group}_marker_genes.tsv", sep="\t", index=False)

    top20 = df.head(20)
    top20.to_csv(results_dir / f"cluster_{group}_top20_marker_genes.tsv", sep="\t", index=False)

    all_marker_tables.append(df)

all_markers = pd.concat(all_marker_tables, axis=0)
all_markers.to_csv(results_dir / "all_cluster_marker_genes.tsv", sep="\t", index=False)

top20_all = all_markers.groupby("cluster", group_keys=False).head(20)
top20_all.to_csv(results_dir / "all_clusters_top20_marker_genes.tsv", sep="\t", index=False)

top10_all = all_markers.groupby("cluster", group_keys=False).head(10)
top10_all.to_csv(results_dir / "all_clusters_top10_marker_genes.tsv", sep="\t", index=False)

# A compact summary: top 5 marker genes for each cluster
summary_rows = []
for group in groups:
    top5 = (
        all_markers.loc[all_markers["cluster"] == group]
        .head(5)["gene"]
        .tolist()
    )
    summary_rows.append({
        "cluster": group,
        "top5_marker_genes": ",".join(top5)
    })

marker_summary = pd.DataFrame(summary_rows)
marker_summary.to_csv(results_dir / "top5_markers_by_cluster.tsv", sep="\t", index=False)

# -----------------------------
# 5. Plot marker genes
# -----------------------------
# Standard Scanpy marker ranking plot
sc.pl.rank_genes_groups(
    adata,
    key="rank_genes_leiden",
    n_genes=10,
    sharey=False,
    show=False
)
plt.savefig(figures_dir / "rank_genes_groups_top10.png", dpi=150, bbox_inches="tight")
plt.close()

# Dotplot of top marker genes
top_marker_genes = []
for group in groups:
    top_marker_genes.extend(
        all_markers.loc[all_markers["cluster"] == group]
        .head(3)["gene"]
        .tolist()
    )

# Remove duplicated genes while preserving order
top_marker_genes_unique = list(dict.fromkeys(top_marker_genes))

sc.pl.dotplot(
    adata,
    var_names=top_marker_genes_unique,
    groupby="leiden",
    use_raw=True,
    show=False
)
plt.savefig(figures_dir / "top3_marker_genes_dotplot.png", dpi=150, bbox_inches="tight")
plt.close()

# UMAP of several classic PBMC marker genes
classic_markers = [
    "IL7R", "CD3D", "CD3E",     # T cells
    "MS4A1", "CD79A",          # B cells
    "LYZ", "S100A8", "S100A9", # Monocytes
    "NKG7", "GNLY",            # NK cells
    "PPBP",                    # Platelets
    "FCER1A", "CST3"           # Dendritic / myeloid related
]

available_classic_markers = [g for g in classic_markers if g in adata.raw.var_names]

with open(results_dir / "available_classic_pbmc_markers.txt", "w") as f:
    for g in available_classic_markers:
        f.write(g + "\n")

if available_classic_markers:
    sc.pl.umap(
        adata,
        color=available_classic_markers,
        use_raw=True,
        show=False
    )
    plt.savefig(figures_dir / "classic_pbmc_marker_umaps.png", dpi=150, bbox_inches="tight")
    plt.close()

# -----------------------------
# 6. Save AnnData
# -----------------------------
output_h5ad = data_dir / "pbmc3k_marker_genes.h5ad"
adata.write(output_h5ad)

# -----------------------------
# 7. Save task summary
# -----------------------------
summary = pd.DataFrame({
    "item": [
        "cells",
        "clusters",
        "groupby",
        "method",
        "use_raw",
        "top_marker_table",
        "classic_markers_available"
    ],
    "value": [
        adata.n_obs,
        len(groups),
        "leiden",
        "wilcoxon",
        True,
        "results/task07/all_clusters_top20_marker_genes.tsv",
        ",".join(available_classic_markers)
    ]
})

summary.to_csv(results_dir / "task07_marker_gene_summary.tsv", sep="\t", index=False)

print("\nTask summary:")
print(summary)

print("\nTop 5 markers by cluster:")
print(marker_summary)

print("\nMain outputs:")
print(f"  {output_h5ad}")
print(f"  {results_dir / 'all_clusters_top20_marker_genes.tsv'}")
print(f"  {results_dir / 'top5_markers_by_cluster.tsv'}")
print(f"  {figures_dir}")
print("Task 07 finished.")
