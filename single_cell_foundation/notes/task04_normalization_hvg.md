# Task 04: Normalization, Log Transformation, and Highly Variable Genes

## Goal

The goal of this task was to normalize the QC-filtered PBMC3k single-cell RNA-seq data, perform log transformation, and identify highly variable genes.

## Input

The input file was:

~~~text
data/pbmc3k_qc_filtered.h5ad
~~~

This file was generated after Task 03 quality control.

## Normalization

Different cells can have different total counts because of differences in RNA capture efficiency and sequencing depth.

Normalization adjusts each cell to a comparable total count scale.

In this task, each cell was normalized to:

~~~text
target_sum = 10000
~~~

The main command was:

~~~python
sc.pp.normalize_total(adata, target_sum=1e4)
~~~

After normalization, the mean and median normalized total counts were both 10000.

## Log Transformation

After normalization, expression values were log-transformed using:

~~~python
sc.pp.log1p(adata)
~~~

This applies:

~~~text
log(x + 1)
~~~

The +1 allows zero values to remain valid.

## Highly Variable Genes

Highly variable genes, or HVGs, are genes that vary strongly across cells after accounting for average expression.

These genes are useful for downstream dimensionality reduction and clustering.

The command was:

~~~python
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)
~~~

## Results

Summary:

~~~text
cells                       2638
genes_after_qc              13714
highly_variable_genes       1838
normalization_target_sum    10000
~~~

A total of 1838 highly variable genes were selected from 13714 genes.

## Data Storage

Before normalization, raw counts were saved in:

~~~text
adata.layers["counts"]
~~~

After normalization but before log transformation, normalized expression values were saved in:

~~~text
adata.layers["normalized"]
~~~

After log transformation:

~~~text
adata.X
~~~

contains log-normalized expression values.

The highly variable gene annotation was saved in:

~~~text
adata.var["highly_variable"]
~~~

## Outputs

Main output object:

- `data/pbmc3k_normalized_hvg.h5ad`

Main result tables:

- `results/task04/task04_normalization_hvg_summary.tsv`
- `results/task04/pbmc3k_hvg_table.tsv`
- `results/task04/pbmc3k_highly_variable_genes.tsv`
- `results/task04/pbmc3k_top30_hvg.tsv`

Main figures:

- `figures/task04/raw_total_counts_hist.png`
- `figures/task04/normalized_total_counts_hist.png`
- `figures/task04/hvg_mean_dispersion_plot.png`

## Important Notes

After this task, `adata.X` is no longer raw counts.

It contains log-normalized expression values.

Raw counts are preserved in:

~~~text
adata.layers["counts"]
~~~

This distinction is important because different downstream steps expect different forms of the expression matrix.

The normalized total counts plot looks very concentrated because all cells were normalized to the same total count scale.

## Next Step

The next task is dimensionality reduction:

- PCA
- nearest-neighbor graph construction
- UMAP visualization
