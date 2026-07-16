# Task 06: Leiden Clustering

## Goal

The goal of this task was to cluster PBMC3k cells using the nearest-neighbor graph generated in Task 05.

## Input

The input file was:

~~~text
data/pbmc3k_pca_neighbors_umap.h5ad
~~~

This file contains:

- PCA coordinates
- nearest-neighbor graph
- UMAP coordinates
- QC metrics
- highly variable gene annotations

## What Leiden Clustering Does

Leiden clustering is a graph-based community detection method.

It does not cluster cells directly from the UMAP plot.

Instead, it uses the nearest-neighbor graph created from PCA space.

Cells that are strongly connected to each other in this graph are assigned to the same cluster.

## Resolution

The resolution parameter controls how coarse or fine the clustering is.

Lower resolution gives fewer and larger clusters.

Higher resolution gives more and smaller clusters.

In this task, four resolutions were tested:

~~~text
0.3
0.5
0.8
1.0
~~~

The final primary clustering used:

~~~text
resolution = 0.8
~~~

The primary cluster labels were stored in:

~~~text
adata.obs["leiden"]
~~~

## Results

The primary clustering result contained 8 clusters:

~~~text
cluster 0    1125 cells
cluster 1     340 cells
cluster 2     222 cells
cluster 3     429 cells
cluster 4     417 cells
cluster 5      56 cells
cluster 6      36 cells
cluster 7      13 cells
~~~

Resolution comparison:

~~~text
resolution 0.3    6 clusters
resolution 0.5    6 clusters
resolution 0.8    8 clusters
resolution 1.0    9 clusters
~~~

Resolution 0.8 was selected as the primary clustering because it gives a more useful level of granularity for PBMC3k marker gene analysis than resolution 0.5, while avoiding excessive fragmentation.

## Important Concept

Cluster is not the same as cell type.

A cluster is an algorithmic group of cells with similar expression patterns.

A cell type is a biological interpretation based on marker genes and prior knowledge.

Therefore, Leiden clusters need to be annotated in later tasks using marker genes.

## Outputs

Main output object:

- `data/pbmc3k_leiden.h5ad`

Main result tables:

- `results/task06/task06_leiden_clustering_summary.tsv`
- `results/task06/leiden_resolution_summary.tsv`
- `results/task06/leiden_primary_cluster_counts.tsv`
- `results/task06/cell_leiden_clusters.tsv`
- `results/task06/leiden_r0_3_cluster_counts.tsv`
- `results/task06/leiden_r0_5_cluster_counts.tsv`
- `results/task06/leiden_r0_8_cluster_counts.tsv`
- `results/task06/leiden_r1_0_cluster_counts.tsv`

Main figures:

- `figures/task06/umap_leiden_primary_r0_8.png`
- `figures/task06/umap_leiden_r0_3.png`
- `figures/task06/umap_leiden_r0_5.png`
- `figures/task06/umap_leiden_r0_8.png`
- `figures/task06/umap_leiden_r1_0.png`
- `figures/task06/leiden_primary_cluster_sizes.png`

## Next Step

The next task is marker gene detection.

Marker genes will be used to interpret the biological cell types represented by each cluster.
