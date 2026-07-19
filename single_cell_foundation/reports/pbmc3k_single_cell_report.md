# PBMC3k Single-cell RNA-seq Analysis Report

## 1. Project Overview

This project analyzed the PBMC3k single-cell RNA-seq dataset using Scanpy.

PBMC means peripheral blood mononuclear cells. The goal of this project was to complete a beginner-level single-cell RNA-seq workflow, from loading the count matrix to annotating major immune cell types.

The workflow included AnnData loading, quality control, normalization, log transformation, highly variable gene selection, PCA, nearest-neighbor graph construction, UMAP visualization, Leiden clustering, marker gene detection, and manual cell type annotation.

## 2. Dataset Summary

| Stage | Cells | Genes |
|---|---:|---:|
| Raw dataset | 2700 | 32738 |
| After quality control | 2638 | 13714 |

The raw PBMC3k dataset contained 2700 cells and 32738 genes. After quality control filtering, 2638 cells and 13714 genes were retained for downstream analysis.

## 3. AnnData Structure

The analysis used the AnnData object as the central data structure.

| AnnData field | Meaning |
|---|---|
| `adata.X` | Main expression matrix |
| `adata.obs` | Cell-level metadata |
| `adata.var` | Gene-level metadata |
| `adata.layers["counts"]` | Raw counts |
| `adata.layers["normalized"]` | Normalized counts before log transformation |
| `adata.raw` | Log-normalized expression before scaling |
| `adata.obsm["X_pca"]` | PCA coordinates |
| `adata.obsm["X_umap"]` | UMAP coordinates |

A key point is that `adata.X` changes during preprocessing. After normalization and log transformation, it no longer represents raw counts. Raw counts were preserved in `adata.layers["counts"]`.

## 4. Quality Control

Quality control was performed to remove low-quality cells and rarely detected genes.

The main cell-level QC metrics were:

| Metric | Meaning |
|---|---|
| `n_genes_by_counts` | Number of detected genes per cell |
| `total_counts` | Total detected UMI/counts per cell |
| `pct_counts_mt` | Percentage of counts from mitochondrial genes |

Human mitochondrial genes were identified by the `MT-` prefix.

The filtering strategy was:

| Filter | Threshold |
|---|---:|
| Minimum genes per cell | 200 |
| Maximum genes per cell | 2500 |
| Maximum mitochondrial percentage | 5% |
| Minimum cells per gene | 3 |

QC results:

| Item | Value |
|---|---:|
| Raw cells | 2700 |
| Filtered cells | 2638 |
| Removed cells | 62 |
| Raw genes | 32738 |
| Filtered genes | 13714 |
| Removed genes | 19024 |
| Mitochondrial genes detected | 13 |

The maximum mitochondrial percentage decreased from about 22.57% before filtering to about 4.99% after filtering. This indicates that cells with high mitochondrial content were successfully removed.

## 5. Normalization and Highly Variable Genes

After quality control, the expression matrix was normalized so that each cell had a comparable total count scale.

Each cell was normalized to a target total count of 10000. The data were then log-transformed using `log(x + 1)`.

Highly variable genes were selected to focus downstream analysis on genes that vary meaningfully across cells.

| Item | Value |
|---|---:|
| Cells after QC | 2638 |
| Genes after QC | 13714 |
| Highly variable genes | 1838 |
| Normalization target sum | 10000 |

After this step, 1838 highly variable genes were selected for dimensionality reduction and clustering.

## 6. PCA, Nearest-Neighbor Graph, and UMAP

The normalized and log-transformed expression data were scaled before PCA.

PCA was used to compress the high-dimensional expression matrix into a smaller set of principal components. These principal components were then used to construct a nearest-neighbor graph.

| Parameter | Value |
|---|---:|
| PCA components computed | 50 |
| PCs used for neighbors | 40 |
| Number of neighbors | 10 |
| UMAP dimensions | 2 |

UMAP was used to visualize cell-cell transcriptomic similarity in two dimensions. Each point in the UMAP represents one cell. Cells close to each other on UMAP usually have similar expression profiles, but UMAP axes and distances should not be over-interpreted as exact biological measurements.

## 7. Leiden Clustering

Leiden clustering was performed on the nearest-neighbor graph.

Several resolutions were tested:

| Resolution | Number of clusters |
|---:|---:|
| 0.3 | 6 |
| 0.5 | 6 |
| 0.8 | 8 |
| 1.0 | 9 |

Resolution 0.8 was selected as the primary clustering result because it provided a useful level of granularity for PBMC marker gene analysis without excessive fragmentation.

Primary cluster sizes:

| Cluster | Number of cells |
|---:|---:|
| 0 | 1125 |
| 1 | 340 |
| 2 | 222 |
| 3 | 429 |
| 4 | 417 |
| 5 | 56 |
| 6 | 36 |
| 7 | 13 |

A cluster is an algorithmic group of cells with similar expression patterns. It is not automatically the same as a biological cell type. Biological interpretation requires marker gene analysis.

## 8. Marker Gene Detection

Marker genes were detected for each Leiden cluster by comparing each cluster against all other cells.

Top marker genes by cluster:

| Cluster | Top marker genes |
|---:|---|
| 0 | LDHB, RPS12, RPS25, RPS3, RPS27 |
| 1 | CD74, CD79A, HLA-DRA, CD79B, HLA-DPB1 |
| 2 | LST1, COTL1, AIF1, FCER1G, FTH1 |
| 3 | NKG7, CST7, GZMA, CTSW, B2M |
| 4 | S100A9, S100A8, LYZ, TYROBP, FTL |
| 5 | MALAT1, RPL32, RPS27, RPL27A, RPS15A |
| 6 | HLA-DPA1, HLA-DPB1, HLA-DRA, HLA-DRB1, CD74 |
| 7 | PF4, SDPR, GNG11, PPBP, NRGN |

Some clusters had very clear marker genes. For example, cluster 1 expressed B cell markers such as `CD79A`, `CD79B`, and `MS4A1`, while cluster 7 expressed platelet markers such as `PF4` and `PPBP`.

Cluster 5 was less clear because its top markers were mainly `MALAT1` and ribosomal genes, which are not specific enough for confident cell type annotation.

## 9. Cell Type Annotation

Manual marker-based annotation was performed using known PBMC marker genes.

Final annotation:

| Cluster | Annotated cell type | Confidence | Evidence genes |
|---:|---|---|---|
| 0 | T cells | Medium | CD3D, CD3E, IL7R, LTB |
| 1 | B cells | High | CD79A, CD79B, MS4A1, CD74, HLA-DRA |
| 2 | Monocytes / myeloid cells | Medium | LST1, COTL1, AIF1, FCER1G, FTH1 |
| 3 | NK / cytotoxic cells | High | NKG7, CST7, GZMA, CTSW, GNLY |
| 4 | CD14+ monocytes | High | S100A8, S100A9, LYZ, TYROBP, FTL |
| 5 | T-like / uncertain | Low | MALAT1, RPL32, RPS27, RPL27A, RPS15A |
| 6 | Dendritic cells | Medium | HLA-DPA1, HLA-DPB1, HLA-DRA, HLA-DRB1, CD74, FCER1A, CST3 |
| 7 | Platelets | High | PF4, PPBP, SDPR, GNG11, NRGN |

Final cell type counts:

| Cell type | Number of cells |
|---|---:|
| T cells | 1125 |
| NK / cytotoxic cells | 429 |
| CD14+ monocytes | 417 |
| B cells | 340 |
| Monocytes / myeloid cells | 222 |
| T-like / uncertain | 56 |
| Dendritic cells | 36 |
| Platelets | 13 |

## 10. Biological Interpretation

The analysis identified several major PBMC populations, including T cells, B cells, monocytes, NK/cytotoxic cells, dendritic cells, and platelets.

The largest population was T cells, followed by NK/cytotoxic cells, CD14+ monocytes, B cells, and another monocyte/myeloid population. Smaller populations included dendritic cells and platelets.

The marker gene patterns were consistent with known PBMC biology. T cell-associated markers such as `CD3D`, `CD3E`, and `IL7R` were enriched in the T cell region. B cell markers such as `MS4A1` and `CD79A` were enriched in the B cell cluster. Monocyte markers such as `LYZ`, `S100A8`, and `S100A9` were enriched in myeloid clusters. Cytotoxic markers such as `NKG7` and `GNLY` supported the NK/cytotoxic annotation. Platelet markers such as `PF4` and `PPBP` supported the platelet annotation.

Cluster 5 was kept as T-like / uncertain because it did not show sufficiently specific marker genes. This is a cautious annotation and should not be over-interpreted.

## 11. Main Output Files

Important result files:

| File | Description |
|---|---|
| `results/task03/pbmc3k_qc_filtering_summary.tsv` | QC filtering summary |
| `results/task04/task04_normalization_hvg_summary.tsv` | Normalization and HVG summary |
| `results/task05/task05_pca_neighbors_umap_summary.tsv` | PCA, neighbors, and UMAP summary |
| `results/task06/leiden_primary_cluster_counts.tsv` | Leiden cluster counts |
| `results/task07/top5_markers_by_cluster.tsv` | Top marker genes by cluster |
| `results/task08/cluster_annotation_summary.tsv` | Final cluster annotation summary |
| `results/task08/cell_type_counts.tsv` | Final cell type counts |

Important figures:

| Figure | Description |
|---|---|
| `figures/task03/raw_pct_counts_mt_hist.png` | Mitochondrial percentage QC |
| `figures/task04/hvg_mean_dispersion_plot.png` | Highly variable gene selection |
| `figures/task05/umap_plain.png` | Basic UMAP visualization |
| `figures/task06/umap_leiden_primary_r0_8.png` | Leiden clustering on UMAP |
| `figures/task07/classic_pbmc_marker_umaps.png` | Classic marker expression on UMAP |
| `figures/task08/umap_cell_type_annotation.png` | Final cell type annotation |
| `figures/task08/annotation_marker_dotplot_by_cluster.png` | Marker gene dotplot by cluster |
| `figures/task08/annotation_marker_dotplot_by_cell_type.png` | Marker gene dotplot by cell type |

## 12. Limitations

This project is an introductory single-cell RNA-seq analysis and has several limitations.

First, the annotation was manual and marker-based. It was not validated using external reference mapping or experimental validation.

Second, cluster 5 remained uncertain because its marker genes were not sufficiently specific.

Third, the analysis used a small tutorial dataset. Real projects often include multiple samples, batches, conditions, donors, and more complex batch correction or integration steps.

Fourth, UMAP was used as a visualization method. UMAP structure should be interpreted together with marker genes, not as direct proof of cell type identity.

## 13. Conclusion

This project successfully completed a basic PBMC3k single-cell RNA-seq workflow using Scanpy.

The workflow progressed from raw AnnData loading to quality control, normalization, highly variable gene selection, PCA, nearest-neighbor graph construction, UMAP visualization, Leiden clustering, marker gene detection, and cell type annotation.

The final marker-based annotation identified major PBMC cell populations, including T cells, B cells, monocytes, NK/cytotoxic cells, dendritic cells, and platelets.

This project demonstrates the core logic of single-cell RNA-seq analysis: start with a count matrix, build a quality-controlled AnnData object, construct a low-dimensional and graph-based representation, cluster cells, identify marker genes, and interpret clusters as biological cell types.
