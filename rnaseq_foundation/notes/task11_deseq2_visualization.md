# Task 11: Visualization of DESeq2 Results

## 1. Task Goal

The goal of this task was to visualize and summarize the DESeq2 differential expression results from Task 10.

The main comparison was:

- dexamethasone-treated samples
- untreated samples

## 2. Input Files

Main input files:

- results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv
- results/deseq2/dds_deseq2.rds
- results/deseq2/normalized_counts.tsv
- metadata/samples.csv

## 3. Volcano Plot

A volcano plot was generated to show the relationship between effect size and statistical significance.

Output:

- results/figures/task11_volcano_plot.pdf

Genes with padj < 0.05 and log2FoldChange >= 1 were labeled as upregulated.

Genes with padj < 0.05 and log2FoldChange <= -1 were labeled as downregulated.

## 4. Candidate Gene Table

Candidate genes were selected using the following criteria:

- padj < 0.05
- absolute log2FoldChange >= 1
- baseMean >= 100

Output files:

- results/task11/candidate_genes_padj_0.05_abs_log2fc_1_basemean_100.tsv
- results/task11/top50_candidate_genes.tsv

These genes are useful for downstream biological interpretation and potential experimental validation.

## 5. Top 50 Heatmap

A heatmap was generated using the top 50 candidate genes.

Output:

- results/figures/task11_top50_candidate_gene_heatmap.pdf

The heatmap shows whether the strongest differentially expressed genes separate dexamethasone-treated samples from untreated samples.

## 6. Selected Gene Expression Plot

Selected genes were visualized across paired cell lines.

Output:

- results/figures/task11_selected_gene_expression.pdf

Selected genes included:

- ZBTB16
- DUSP1
- PER1
- VCAM1
- SAMHD1
- GPX3
- SERPINA3
- ADAMTS1

This plot helps check whether important genes show consistent treatment-related changes across cell lines.

## 7. Initial Interpretation

The visualization results provide a clearer view of the dexamethasone response.

Strongly upregulated genes include known or plausible glucocorticoid-responsive genes such as ZBTB16, DUSP1, and PER1.

VCAM1 showed strong downregulation after dexamethasone treatment.

## 8. Next Step

The next step is Task 12: pathway enrichment analysis.

This will include GO enrichment and/or GSEA to identify biological processes affected by dexamethasone treatment.
