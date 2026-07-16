# Task 01: Project Setup

## Goal

Initialize the single-cell RNA-seq foundation project.

## Project Name

single_cell_foundation

## Motivation

This project follows the completed bulk RNA-seq foundation project.

The goal is to learn how single-cell RNA-seq analysis differs from bulk RNA-seq.

Bulk RNA-seq summarizes expression at the sample level:

~~~text
gene × sample
~~~

Single-cell RNA-seq measures expression at the cell level:

~~~text
gene × cell
~~~

## Planned Dataset

PBMC 3k from 10x Genomics.

PBMC means peripheral blood mononuclear cells.

This dataset is suitable for beginner-level scRNA-seq analysis because:

- it is small
- it is well documented
- it is widely used in Scanpy and Seurat tutorials
- it contains recognizable immune cell types

## Planned Workflow

~~~text
Count matrix
↓
QC
↓
Normalization
↓
Highly variable genes
↓
PCA
↓
Neighbors graph
↓
UMAP
↓
Clustering
↓
Marker genes
↓
Cell type annotation
↓
Report
~~~

## Notes

This project will use Python and Scanpy instead of R/DESeq2.

The goal is not only to run functions, but to understand what each object and plot means.
