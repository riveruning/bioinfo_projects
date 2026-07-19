# Task 08: Cell Type Annotation

## Goal

The goal of this task was to annotate PBMC3k Leiden clusters with biological cell type labels based on marker genes.

## Input

The input file was:

~~~text
data/pbmc3k_marker_genes.h5ad
~~~

This file contains:

- QC-filtered cells
- normalized and log-transformed expression values
- PCA coordinates
- nearest-neighbor graph
- UMAP coordinates
- Leiden cluster labels
- marker gene ranking results

## Annotation Strategy

Cell type annotation was performed manually using known PBMC marker genes.

The logic was:

~~~text
Leiden cluster
↓
cluster marker genes
↓
known immune cell marker genes
↓
cell type label
~~~

This is a marker-based manual annotation, not an automatic ground truth label.

## Cluster Annotation

The first-pass annotation was:

~~~text
cluster 0    T cells
cluster 1    B cells
cluster 2    Monocytes / myeloid cells
cluster 3    NK / cytotoxic cells
cluster 4    CD14+ monocytes
cluster 5    T-like / uncertain
cluster 6    Dendritic cells
cluster 7    Platelets
~~~

## Evidence

Examples of marker evidence:

- T cells: `CD3D`, `CD3E`, `IL7R`, `LTB`
- B cells: `MS4A1`, `CD79A`, `CD79B`, `CD74`
- Monocytes: `LYZ`, `LST1`, `COTL1`, `AIF1`, `S100A8`, `S100A9`
- NK / cytotoxic cells: `NKG7`, `GNLY`, `GZMA`, `CST7`
- Dendritic cells: `FCER1A`, `CST3`, `HLA-DRA`, `HLA-DPA1`, `HLA-DPB1`
- Platelets: `PF4`, `PPBP`, `GNG11`, `SDPR`

## Important Caveat

Cluster 5 was annotated as:

~~~text
T-like / uncertain
~~~

because its top marker genes were mainly `MALAT1` and ribosomal genes such as `RPL` and `RPS` genes.

These genes are not sufficiently cell-type-specific.

Therefore, cluster 5 should be treated cautiously.

## Outputs

Main output object:

- `data/pbmc3k_annotated.h5ad`

Main result tables:

- `results/task08/task08_cell_type_annotation_summary.tsv`
- `results/task08/cluster_cell_type_annotation.tsv`
- `results/task08/cluster_annotation_summary.tsv`
- `results/task08/cell_type_counts.tsv`
- `results/task08/cell_cluster_cell_type.tsv`
- `results/task08/available_annotation_marker_panel.txt`

Main figures:

- `figures/task08/umap_leiden.png`
- `figures/task08/umap_cell_type_annotation.png`
- `figures/task08/annotation_marker_dotplot_by_cluster.png`
- `figures/task08/annotation_marker_dotplot_by_cell_type.png`
- `figures/task08/selected_marker_umaps.png`

## Interpretation

The annotated UMAP connects computational clustering with biological interpretation.

The clusters are no longer only numerical labels.

They now have first-pass biological labels supported by marker gene expression.

## Next Step

The next task is to write a final project report summarizing the full single-cell RNA-seq workflow and biological interpretation.
