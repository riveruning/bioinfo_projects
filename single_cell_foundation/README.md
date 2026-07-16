# Single-cell RNA-seq Foundation Project

This project is a beginner-friendly single-cell RNA-seq analysis project using Python and Scanpy.

The goal is to learn the standard scRNA-seq workflow from a count matrix to quality control, normalization, dimensionality reduction, clustering, marker gene detection, cell type annotation, and biological interpretation.

## Project Goal

The main question is:

> How can we identify major cell populations from a single-cell RNA-seq count matrix?

## Planned Workflow

~~~text
Count matrix
↓
AnnData object
↓
Quality control
↓
Normalization and log transform
↓
Highly variable genes
↓
PCA
↓
Neighbors graph
↓
UMAP
↓
Leiden clustering
↓
Marker gene detection
↓
Cell type annotation
↓
Biological interpretation
~~~

## Dataset

The first dataset will be PBMC 3k from 10x Genomics.

PBMC means peripheral blood mononuclear cells.

This dataset is commonly used for introductory scRNA-seq tutorials because it is small, well documented, and contains recognizable immune cell populations.

## Main Tools

- Python
- Scanpy
- AnnData
- pandas
- matplotlib

## Directory Structure

~~~text
data/       raw and processed data files
metadata/   sample or cell metadata
scripts/    analysis scripts
results/    result tables and processed objects
figures/    plots and visual outputs
notes/      task-by-task learning notes
reports/    final reports
~~~

## Current Status

Project initialized.

Next task:

- Set up the Scanpy environment
- Load the PBMC 3k dataset
- Inspect the AnnData object
