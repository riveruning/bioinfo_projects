# Task 07: Marker Gene Detection

## Goal

The goal of this task was to identify marker genes for each Leiden cluster in the PBMC3k single-cell RNA-seq dataset.

## Input

The input file was:

~~~text
data/pbmc3k_leiden.h5ad
~~~

This file contains:

- QC-filtered cells
- normalized and log-transformed expression values
- PCA coordinates
- nearest-neighbor graph
- UMAP coordinates
- Leiden cluster labels

## What Marker Genes Are

Marker genes are genes that are relatively highly expressed in one cluster compared with other clusters.

They help connect algorithmic clusters to biological cell identities.

For example, in PBMC data:

- `CD3D`, `CD3E`, and `IL7R` are often associated with T cells
- `MS4A1` and `CD79A` are often associated with B cells
- `LYZ`, `S100A8`, and `S100A9` are often associated with monocytes
- `NKG7` and `GNLY` are often associated with NK cells
- `PPBP` is often associated with platelets

## Method

Marker genes were detected using:

~~~python
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    use_raw=True,
    key_added="rank_genes_leiden"
)
~~~

The comparison was performed for each Leiden cluster against all other cells.

## Why use_raw=True Was Used

In Task 05, expression values in `adata.X` were scaled for PCA.

Scaled expression values are useful for dimensionality reduction but are less suitable for direct biological interpretation.

The object `adata.raw` stores log-normalized expression values before scaling.

Therefore, marker gene detection used:

~~~text
use_raw = True
~~~

## Outputs

Main output object:

- `data/pbmc3k_marker_genes.h5ad`

Main result tables:

- `results/task07/task07_marker_gene_summary.tsv`
- `results/task07/top5_markers_by_cluster.tsv`
- `results/task07/all_clusters_top10_marker_genes.tsv`
- `results/task07/all_clusters_top20_marker_genes.tsv`
- `results/task07/all_cluster_marker_genes.tsv`
- `results/task07/cluster_*_marker_genes.tsv`
- `results/task07/cluster_*_top20_marker_genes.tsv`
- `results/task07/available_classic_pbmc_markers.txt`

Main figures:

- `figures/task07/rank_genes_groups_top10.png`
- `figures/task07/top3_marker_genes_dotplot.png`
- `figures/task07/classic_pbmc_marker_umaps.png`

## Important Concept

Cluster is not the same as cell type.

Marker genes provide evidence for interpreting each cluster as a possible cell type.

Formal cell type annotation will be performed in the next task.

## Next Step

The next task is cell type annotation.
