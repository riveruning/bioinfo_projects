# Task 05: PCA, Nearest-Neighbor Graph, and UMAP

## Goal

The goal of this task was to perform dimensionality reduction and construct a cell-cell similarity graph for the QC-filtered and normalized PBMC3k single-cell RNA-seq dataset.

## Input

The input file was:

~~~text
data/pbmc3k_normalized_hvg.h5ad
~~~

This file was generated in Task 04.

It contains:

- QC-filtered cells
- normalized and log-transformed expression values
- highly variable gene annotations

## Scaling

Before PCA, expression values were scaled using:

~~~python
sc.pp.scale(adata, max_value=10)
~~~

Scaling makes genes more comparable for PCA.

The parameter `max_value=10` clips extreme scaled values to reduce the effect of outliers.

## PCA

PCA means principal component analysis.

It compresses high-dimensional gene expression information into a smaller number of principal components.

In this task, PCA was computed using 50 components and highly variable genes.

## Nearest-Neighbor Graph

A nearest-neighbor graph connects each cell to cells with similar PCA representations.

The graph was computed using:

~~~text
n_neighbors = 10
n_pcs = 40
~~~

This graph is the basis for UMAP visualization and later Leiden clustering.

## UMAP

UMAP is a nonlinear dimensionality reduction method used for visualization.

It projects cells into two dimensions while trying to preserve local neighborhood relationships.

## Results

Main settings:

~~~text
cells                  2638
genes                  13714
highly_variable_genes  1838
pca_n_comps            50
neighbors_n_neighbors  10
neighbors_n_pcs        40
umap_dimensions        2
~~~

The UMAP showed several separated cell groups, suggesting that the PBMC3k dataset contains multiple transcriptionally distinct cell populations.

The mitochondrial percentage did not show a strong abnormal pattern across the UMAP, indicating that the previous QC filtering was acceptable.

## Outputs

Main output object:

- `data/pbmc3k_pca_neighbors_umap.h5ad`

Main result tables:

- `results/task05/task05_pca_neighbors_umap_summary.tsv`
- `results/task05/pca_variance_ratio.tsv`
- `results/task05/pca_coordinates_first10.tsv`
- `results/task05/umap_coordinates.tsv`

Main figures:

- `figures/task05/pca_variance_ratio.png`
- `figures/task05/umap_plain.png`
- `figures/task05/umap_n_genes_by_counts.png`
- `figures/task05/umap_total_counts.png`
- `figures/task05/umap_pct_counts_mt.png`

## Important Notes

PCA is used as an intermediate representation for downstream analysis.

UMAP is mainly a visualization tool.

Distances on UMAP should not be over-interpreted as exact biological distances.

The nearest-neighbor graph created in this task will be used for Leiden clustering in the next task.

## Next Step

The next task is Leiden clustering.
