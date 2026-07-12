# Task 10: DESeq2 Differential Expression Analysis

## 1. Task Goal

The goal of this task was to perform differential expression analysis for the airway RNA-seq dataset.

The main comparison was:

- dexamethasone-treated samples
- untreated samples

## 2. Input Files

The main input was the tximport gene-level object:

- results/tximport/txi_salmon_gene_level.rds

Sample metadata was read from:

- metadata/samples.csv

Gene annotation was read from:

- results/tximport/gene_annotation.tsv

## 3. Experimental Design

The airway dataset contains four cell lines.

Each cell line has one untreated sample and one dexamethasone-treated sample.

Therefore, the DESeq2 design used was:

design = ~ cell_line + dex

This model accounts for baseline differences between cell lines and then estimates the dexamethasone treatment effect.

## 4. Filtering

Low-count genes were filtered before running DESeq2.

Genes were kept if they had count >= 10 in at least 2 samples.

This reduced the number of genes used for statistical testing.

## 5. Summary

The DESeq2 summary was:

- samples: 8
- genes before filtering: 82323
- genes after filtering: 20407
- significant genes with padj < 0.05: 3711
- upregulated genes with padj < 0.05: 2077
- downregulated genes with padj < 0.05: 1634
- significant genes with padj < 0.05 and |log2FC| >= 1: 988

## 6. Outputs

Main DESeq2 output files:

- results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv
- results/deseq2/significant_genes_padj_0.05.tsv
- results/deseq2/significant_genes_padj_0.05_abs_log2fc_1.tsv
- results/deseq2/normalized_counts.tsv
- results/deseq2/task10_summary.tsv

Main figures:

- results/figures/task10_pca_plot.pdf
- results/figures/task10_ma_plot.pdf
- results/figures/task10_sample_distance_heatmap.pdf

## 7. Annotation Issue and Fix

An initial version of the result table contained duplicated gene rows after merging gene annotation.

This happened because the transcript-derived gene annotation file contained multiple gene_type values for the same gene_id.

The script was fixed by collapsing annotation to one row per gene_id before merging with DESeq2 results.

After this fix, each gene_id appears only once in the DESeq2 result table.

## 8. Initial Interpretation

Positive log2FoldChange means the gene is upregulated after dexamethasone treatment.

Negative log2FoldChange means the gene is downregulated after dexamethasone treatment.

Top-ranked genes included NEXN, ZBTB16, DUSP1, SORT1, VCAM1, SAMHD1, MT2A, ADAMTS1, SERPINA3, MAOA, and PER1.

Several of these genes show strong responses to dexamethasone treatment.

## 9. Next Step

The next step is Task 11: visualize and interpret the differential expression results.

This will include volcano plots, top gene heatmaps, and checking known dexamethasone-responsive genes.
