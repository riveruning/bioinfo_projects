import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

print("Task 08: cell type annotation")

# -----------------------------
# 1. Paths
# -----------------------------
input_h5ad = Path("data/pbmc3k_marker_genes.h5ad")
data_dir = Path("data")
results_dir = Path("results/task08")
figures_dir = Path("figures/task08")

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
    raise ValueError("Missing adata.raw. Marker visualization needs log-normalized raw data.")

# -----------------------------
# 3. Manual annotation based on Task 07 marker genes
# -----------------------------
cluster_to_cell_type = {
    "0": "T cells",
    "1": "B cells",
    "2": "Monocytes / myeloid cells",
    "3": "NK / cytotoxic cells",
    "4": "CD14+ monocytes",
    "5": "T-like / uncertain",
    "6": "Dendritic cells",
    "7": "Platelets",
}

cluster_evidence = [
    {
        "cluster": "0",
        "cell_type": "T cells",
        "evidence_genes": "CD3D,CD3E,IL7R,LTB",
        "confidence": "medium",
        "comment": "T cell markers are present, although ribosomal genes also rank highly."
    },
    {
        "cluster": "1",
        "cell_type": "B cells",
        "evidence_genes": "CD79A,CD79B,MS4A1,CD74,HLA-DRA",
        "confidence": "high",
        "comment": "B cell receptor and MHC-II genes support B cell identity."
    },
    {
        "cluster": "2",
        "cell_type": "Monocytes / myeloid cells",
        "evidence_genes": "LST1,COTL1,AIF1,FCER1G,FTH1",
        "confidence": "medium",
        "comment": "Myeloid and monocyte-associated genes are enriched."
    },
    {
        "cluster": "3",
        "cell_type": "NK / cytotoxic cells",
        "evidence_genes": "NKG7,CST7,GZMA,CTSW,GNLY",
        "confidence": "high",
        "comment": "Cytotoxic markers support NK or cytotoxic lymphocyte identity."
    },
    {
        "cluster": "4",
        "cell_type": "CD14+ monocytes",
        "evidence_genes": "S100A8,S100A9,LYZ,TYROBP,FTL",
        "confidence": "high",
        "comment": "Classic inflammatory monocyte markers are enriched."
    },
    {
        "cluster": "5",
        "cell_type": "T-like / uncertain",
        "evidence_genes": "MALAT1,RPL32,RPS27,RPL27A,RPS15A",
        "confidence": "low",
        "comment": "Top markers are not cell-type specific. This cluster should be treated cautiously."
    },
    {
        "cluster": "6",
        "cell_type": "Dendritic cells",
        "evidence_genes": "HLA-DPA1,HLA-DPB1,HLA-DRA,HLA-DRB1,CD74,FCER1A,CST3",
        "confidence": "medium",
        "comment": "MHC-II and dendritic/myeloid-related markers support dendritic cell identity."
    },
    {
        "cluster": "7",
        "cell_type": "Platelets",
        "evidence_genes": "PF4,PPBP,SDPR,GNG11,NRGN",
        "confidence": "high",
        "comment": "Platelet markers are strongly enriched."
    },
]

annotation_df = pd.DataFrame(cluster_evidence)
annotation_df.to_csv(results_dir / "cluster_cell_type_annotation.tsv", sep="\t", index=False)

# Add cell type labels to adata.obs
adata.obs["leiden_str"] = adata.obs["leiden"].astype(str)
adata.obs["cell_type"] = adata.obs["leiden_str"].map(cluster_to_cell_type)
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

# -----------------------------
# 4. Save count tables
# -----------------------------
cluster_counts = adata.obs["leiden"].value_counts().sort_index()
cluster_counts_df = cluster_counts.rename_axis("cluster").reset_index(name="n_cells")
cluster_counts_df.to_csv(results_dir / "cluster_counts.tsv", sep="\t", index=False)

cell_type_counts = adata.obs["cell_type"].value_counts()
cell_type_counts_df = cell_type_counts.rename_axis("cell_type").reset_index(name="n_cells")
cell_type_counts_df.to_csv(results_dir / "cell_type_counts.tsv", sep="\t", index=False)

cell_annotation_df = adata.obs[["leiden", "cell_type"]].copy()
cell_annotation_df.to_csv(results_dir / "cell_cluster_cell_type.tsv", sep="\t")

# Merge cluster counts with annotation evidence
cluster_annotation_summary = annotation_df.merge(
    cluster_counts_df,
    on="cluster",
    how="left"
)

cluster_annotation_summary.to_csv(
    results_dir / "cluster_annotation_summary.tsv",
    sep="\t",
    index=False
)

# -----------------------------
# 5. Marker panel for checking annotation
# -----------------------------
marker_panel = {
    "T cells": ["CD3D", "CD3E", "IL7R", "LTB", "CCR7"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD74"],
    "Monocytes": ["LYZ", "LST1", "COTL1", "AIF1", "S100A8", "S100A9", "CD14", "FCGR3A", "FCER1G"],
    "NK / cytotoxic": ["NKG7", "GNLY", "GZMA", "GZMB", "PRF1", "CST7"],
    "Dendritic": ["FCER1A", "CST3", "HLA-DRA", "HLA-DPA1", "HLA-DPB1"],
    "Platelets": ["PF4", "PPBP", "GNG11", "SDPR", "NRGN"],
}

available_marker_panel = {}
for group, genes in marker_panel.items():
    available_genes = [g for g in genes if g in adata.raw.var_names]
    if available_genes:
        available_marker_panel[group] = available_genes

with open(results_dir / "available_annotation_marker_panel.txt", "w") as f:
    for group, genes in available_marker_panel.items():
        f.write(group + "\t" + ",".join(genes) + "\n")

# -----------------------------
# 6. Figures
# -----------------------------

# UMAP by Leiden cluster
sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    title="Leiden clusters",
    show=False
)
plt.savefig(figures_dir / "umap_leiden.png", dpi=150, bbox_inches="tight")
plt.close()

# UMAP by annotated cell type
sc.pl.umap(
    adata,
    color="cell_type",
    legend_loc="right margin",
    title="Annotated cell types",
    show=False
)
plt.savefig(figures_dir / "umap_cell_type_annotation.png", dpi=150, bbox_inches="tight")
plt.close()

# Dotplot by Leiden cluster
sc.pl.dotplot(
    adata,
    var_names=available_marker_panel,
    groupby="leiden",
    use_raw=True,
    show=False
)
plt.savefig(figures_dir / "annotation_marker_dotplot_by_cluster.png", dpi=150, bbox_inches="tight")
plt.close()

# Dotplot by annotated cell type
sc.pl.dotplot(
    adata,
    var_names=available_marker_panel,
    groupby="cell_type",
    use_raw=True,
    show=False
)
plt.savefig(figures_dir / "annotation_marker_dotplot_by_cell_type.png", dpi=150, bbox_inches="tight")
plt.close()

# UMAPs for selected key marker genes
selected_markers = [
    "CD3D", "CD3E", "IL7R",
    "MS4A1", "CD79A",
    "LYZ", "S100A8", "S100A9",
    "NKG7", "GNLY",
    "FCER1A", "CST3",
    "PF4", "PPBP"
]

available_selected_markers = [g for g in selected_markers if g in adata.raw.var_names]

if available_selected_markers:
    sc.pl.umap(
        adata,
        color=available_selected_markers,
        use_raw=True,
        show=False
    )
    plt.savefig(figures_dir / "selected_marker_umaps.png", dpi=150, bbox_inches="tight")
    plt.close()

# -----------------------------
# 7. Save annotated AnnData
# -----------------------------
output_h5ad = data_dir / "pbmc3k_annotated.h5ad"
adata.write(output_h5ad)

# -----------------------------
# 8. Save task summary
# -----------------------------
summary = pd.DataFrame({
    "item": [
        "cells",
        "clusters",
        "annotated_cell_types",
        "annotation_type",
        "uncertain_clusters",
        "output_h5ad"
    ],
    "value": [
        adata.n_obs,
        adata.obs["leiden"].nunique(),
        adata.obs["cell_type"].nunique(),
        "manual_marker_based_annotation",
        "5",
        str(output_h5ad)
    ]
})

summary.to_csv(results_dir / "task08_cell_type_annotation_summary.tsv", sep="\t", index=False)

print("\nTask summary:")
print(summary)

print("\nCluster annotation summary:")
print(cluster_annotation_summary)

print("\nCell type counts:")
print(cell_type_counts_df)

print("\nMain outputs:")
print(f"  {output_h5ad}")
print(f"  {results_dir / 'cluster_cell_type_annotation.tsv'}")
print(f"  {results_dir / 'cell_type_counts.tsv'}")
print(f"  {figures_dir}")
print("Task 08 finished.")
