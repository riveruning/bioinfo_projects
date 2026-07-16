# Task 02: Load PBMC3k and Inspect AnnData

## Goal

The goal of this task was to set up the Scanpy workflow and load the PBMC3k single-cell RNA-seq dataset.

## Dataset

The dataset used in this task was PBMC3k.

PBMC means peripheral blood mononuclear cells.

This is a small 10x Genomics single-cell RNA-seq dataset commonly used for beginner-level single-cell RNA-seq tutorials.

## Main Concepts

The main object used by Scanpy is AnnData.

AnnData stores single-cell data in several connected tables:

- `adata.X`: expression matrix
- `adata.obs`: cell-level metadata
- `adata.var`: gene-level metadata

## Matrix Shape

The PBMC3k AnnData object had the shape:

~~~text
2700 cells × 32738 genes
~~~

In Scanpy, the expression matrix is usually stored as:

~~~text
cells × genes
~~~

This is different from bulk RNA-seq, where expression matrices are often discussed as:

~~~text
genes × samples
~~~

## adata.X

`adata.X` stores the expression matrix.

The matrix type was:

~~~text
scipy.sparse._csr.csr_matrix
~~~

This means the expression matrix is stored as a sparse matrix.

Sparse matrix storage is useful because most gene-by-cell entries in single-cell RNA-seq are zero.

## adata.obs

`adata.obs` stores cell-level information.

At this raw data stage, `adata.obs` did not contain extra metadata columns.

The row index of `adata.obs` contains cell barcodes, such as:

~~~text
AAACATACAACCAC-1
AAACATTGAGCTAC-1
AAACATTGATCAGC-1
~~~

These are 10x cell barcodes.

They are the names or IDs of individual cells in the dataset.

## adata.var

`adata.var` stores gene-level information.

The row index contains gene symbols.

The column `gene_ids` contains Ensembl gene IDs.

Examples:

~~~text
MIR1302-10    ENSG00000243485
FAM138A       ENSG00000237613
OR4F5         ENSG00000186092
~~~

## Important Terminology

- 10x Genomics: a company and technology platform for single-cell sequencing
- cell barcode: a DNA barcode used to label which cell an RNA molecule came from
- index: the row name of a table
- `adata.obs` index: cell barcodes
- `adata.var` index: gene symbols
- h5ad: file format used to save an AnnData object

## Outputs

The script generated:

- `data/pbmc3k_raw.h5ad`
- `results/task02/pbmc3k_basic_summary.tsv`
- `results/task02/pbmc3k_obs_head.tsv`
- `results/task02/pbmc3k_var_head.tsv`
- `results/task02_load_pbmc3k.log`

## Interpretation

At this stage, no filtering, normalization, clustering, or cell type annotation has been performed.

The purpose was only to load the raw single-cell count matrix and understand the AnnData structure.

## Next Step

The next task is single-cell quality control.

This will include calculating:

- number of genes detected per cell
- total UMI counts per cell
- mitochondrial gene percentage
