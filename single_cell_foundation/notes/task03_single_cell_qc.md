# Task 03: Single-cell Quality Control

## Goal

The goal of this task was to calculate quality control metrics for the PBMC3k single-cell RNA-seq dataset and filter low-quality cells and rarely detected genes.

## Input

The input file was:

~~~text
data/pbmc3k_raw.h5ad
~~~

This file contains the raw PBMC3k AnnData object generated in Task 02.

## Main QC Metrics

### n_genes_by_counts

`n_genes_by_counts` is the number of genes detected in each cell.

Cells with very few detected genes may represent low-quality cells, empty droplets, or failed capture events.

Cells with extremely many detected genes may represent doublets.

### total_counts

`total_counts` is the total UMI or count number detected in each cell.

It roughly reflects the total amount of captured RNA molecules for a cell.

### pct_counts_mt

`pct_counts_mt` is the percentage of counts coming from mitochondrial genes.

Human mitochondrial genes usually start with:

~~~text
MT-
~~~

High mitochondrial percentage can indicate damaged or stressed cells.

## Filtering Strategy

Cells were kept if they satisfied:

~~~text
n_genes_by_counts >= 200
n_genes_by_counts < 2500
pct_counts_mt < 5
~~~

Genes were kept if they were detected in at least 3 cells:

~~~text
min_cells = 3
~~~

## Results

Filtering summary:

~~~text
raw_cells        2700
raw_genes        32738
filtered_cells   2638
filtered_genes   13714
removed_cells    62
removed_genes    19024
mitochondrial_genes_detected    13
~~~

The filtering removed a small number of low-quality or abnormal cells and many rarely detected genes.

This is expected for single-cell RNA-seq data.

## Interpretation

The maximum mitochondrial percentage decreased from about 22.57% before filtering to about 4.99% after filtering.

The maximum number of detected genes per cell decreased from 3422 to 2455.

These changes indicate that cells with high mitochondrial percentage or unusually high detected gene counts were removed.

The filtered dataset contains:

~~~text
2638 cells × 13714 genes
~~~

This filtered AnnData object will be used for normalization and downstream analysis.

## Outputs

Main result files:

- `results/task03/pbmc3k_qc_filtering_summary.tsv`
- `results/task03/pbmc3k_cell_qc_metrics_raw.tsv`
- `results/task03/pbmc3k_cell_qc_metrics_filtered.tsv`
- `results/task03/pbmc3k_qc_summary_raw.tsv`
- `results/task03/pbmc3k_qc_summary_filtered.tsv`

Main figures:

- `figures/task03/raw_n_genes_by_counts_hist.png`
- `figures/task03/raw_total_counts_hist.png`
- `figures/task03/raw_pct_counts_mt_hist.png`
- `figures/task03/raw_total_counts_vs_pct_counts_mt.png`
- `figures/task03/raw_total_counts_vs_n_genes.png`
- `figures/task03/filtered_n_genes_by_counts_hist.png`
- `figures/task03/filtered_total_counts_hist.png`
- `figures/task03/filtered_pct_counts_mt_hist.png`

Filtered AnnData object:

- `data/pbmc3k_qc_filtered.h5ad`

## Important Notes

QC thresholds are not universal.

The thresholds used here are suitable for the PBMC3k beginner dataset, but real projects require threshold selection based on the data distribution, tissue type, platform, and biological question.

## Next Step

The next task is normalization, log transformation, and highly variable gene selection.
