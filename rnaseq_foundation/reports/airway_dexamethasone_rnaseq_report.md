# Airway RNA-seq Analysis: Dexamethasone Response in Airway Smooth Muscle Cells

## 1. Project Overview

This project analyzed a public bulk RNA-seq dataset of human airway smooth muscle cells treated with dexamethasone.

The goal was to learn and reproduce a complete RNA-seq analysis workflow, starting from raw FASTQ files and ending with differential expression analysis, visualization, enrichment analysis, and biological interpretation.

The main biological question was:

> Which genes and biological processes are affected by dexamethasone treatment in airway smooth muscle cells?

## 2. Dataset

The dataset contains paired-end bulk RNA-seq samples from human airway smooth muscle cell lines.

Eight samples were analyzed:

| Sample ID | Cell line | Treatment |
|---|---|---|
| SRR1039508 | N61311 | untreated |
| SRR1039509 | N61311 | dexamethasone treated |
| SRR1039512 | N052611 | untreated |
| SRR1039513 | N052611 | dexamethasone treated |
| SRR1039516 | N080611 | untreated |
| SRR1039517 | N080611 | dexamethasone treated |
| SRR1039520 | N061011 | untreated |
| SRR1039521 | N061011 | dexamethasone treated |

The experimental design includes four cell lines, each with one untreated and one dexamethasone-treated sample.

## 3. Analysis Workflow

The analysis workflow was:

~~~text
Raw FASTQ files
↓
FastQC / MultiQC quality control
↓
Salmon transcript-level quantification
↓
tximport transcript-to-gene summarization
↓
DESeq2 differential expression analysis
↓
Volcano plot, heatmap, selected gene expression plots
↓
GO enrichment and GSEA
↓
Biological interpretation and candidate gene selection
~~~

## 4. Software and Tools

Major tools used:

- FastQC: raw sequencing quality control
- MultiQC: multi-sample QC report aggregation
- Salmon: transcript-level RNA-seq quantification
- tximport: transcript-to-gene summarization
- DESeq2: differential expression analysis
- clusterProfiler: GO enrichment and GSEA
- org.Hs.eg.db: human gene ID annotation
- ggplot2 / pheatmap: visualization

## 5. Quality Control

FastQC and MultiQC were used to assess raw FASTQ quality across all samples.

The full MultiQC report is available at:

- `results/qc/multiqc_all_samples/multiqc_report.html`

Overall, the sequencing data were suitable for downstream analysis.

## 6. Salmon Quantification

Salmon was used to quantify transcript abundance against the GENCODE transcriptome reference.

All eight samples showed high mapping rates:

| Sample | Mapping rate |
|---|---:|
| SRR1039508 | 95.42% |
| SRR1039509 | 94.70% |
| SRR1039512 | 96.10% |
| SRR1039513 | 96.54% |
| SRR1039516 | 95.62% |
| SRR1039517 | 96.02% |
| SRR1039520 | 96.03% |
| SRR1039521 | 96.28% |

These high mapping rates suggest that the FASTQ files, reference transcriptome, and quantification setup were appropriate.

Main output:

- `results/salmon_summary/all_samples_mapping_rate.tsv`

## 7. tximport Gene-level Matrix

Salmon transcript-level quantification files were imported with tximport and summarized to gene-level matrices.

Summary:

| Item | Value |
|---|---:|
| Samples | 8 |
| Genes in count matrix | 82323 |
| Transcripts in quant.sf | 654828 |
| Transcripts in tx2gene mapping | 654828 |

Main outputs:

- `results/tximport/gene_counts.tsv`
- `results/tximport/gene_tpm.tsv`
- `results/tximport/txi_salmon_gene_level.rds`
- `results/tximport/tximport_summary.tsv`

## 8. DESeq2 Differential Expression Analysis

DESeq2 was used to compare dexamethasone-treated samples against untreated samples.

Because the dataset includes four paired cell lines, the model controlled for cell line differences:

~~~r
design = ~ cell_line + dex
~~~

This design estimates the dexamethasone treatment effect after accounting for baseline differences between cell lines.

DESeq2 summary:

| Item | Value |
|---|---:|
| Samples | 8 |
| Genes before filtering | 82323 |
| Genes after filtering | 20407 |
| Significant genes, padj < 0.05 | 3711 |
| Upregulated significant genes | 2077 |
| Downregulated significant genes | 1634 |
| Significant genes with padj < 0.05 and abs(log2FC) >= 1 | 988 |

Main output:

- `results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv`

## 9. Differential Expression Visualization

The DESeq2 results were visualized using:

- PCA plot
- MA plot
- sample distance heatmap
- volcano plot
- top 50 candidate gene heatmap
- selected gene expression plot

Key figures:

- `results/figures/task10_pca_plot.pdf`
- `results/figures/task10_ma_plot.pdf`
- `results/figures/task10_sample_distance_heatmap.pdf`
- `results/figures/task11_volcano_plot.pdf`
- `results/figures/task11_top50_candidate_gene_heatmap.pdf`
- `results/figures/task11_selected_gene_expression.pdf`

## 10. Key Differentially Expressed Genes

Several genes showed strong dexamethasone-associated expression changes.

Important upregulated genes included:

| Gene | Direction | Approximate interpretation |
|---|---|---|
| ZBTB16 | upregulated | strong glucocorticoid-responsive candidate |
| DUSP1 | upregulated | stress/MAPK-related response gene |
| PER1 | upregulated | hormone/circadian-response-related gene |
| CRISPLD2 | upregulated | known airway glucocorticoid-responsive candidate |
| ADAMTS1 | upregulated | extracellular matrix remodeling candidate |
| GPX3 | upregulated | oxidative stress / antioxidant-response candidate |
| SERPINA3 | upregulated | secreted / inflammatory-response candidate |

Important downregulated genes included:

| Gene | Direction | Approximate interpretation |
|---|---|---|
| VCAM1 | downregulated | adhesion / inflammatory marker |

Selected approximate fold changes:

| Gene | log2FoldChange | Approximate fold change |
|---|---:|---:|
| ZBTB16 | 6.10 | about 68-fold upregulated |
| DUSP1 | 2.97 | about 7.8-fold upregulated |
| PER1 | 3.03 | about 8.2-fold upregulated |
| CRISPLD2 | 2.61 | about 6.1-fold upregulated |
| VCAM1 | -3.59 | about 12-fold downregulated |

For downregulated genes, fold change should be described as lower expression rather than negative expression.

## 11. GO Enrichment Analysis

GO Biological Process enrichment analysis was performed separately for upregulated and downregulated genes.

Input criteria:

- padj < 0.05
- abs(log2FoldChange) >= 1

Summary:

| Item | Value |
|---|---:|
| Upregulated Entrez genes used for ORA | 485 |
| Downregulated Entrez genes used for ORA | 433 |
| GO BP terms enriched in upregulated genes | 77 |
| GO BP terms enriched in downregulated genes | 418 |

Upregulated genes were enriched for processes including:

- muscle system process
- cell-substrate adhesion
- muscle adaptation
- angiogenesis
- cellular response to glucocorticoid stimulus
- reactive oxygen species metabolic process
- extracellular matrix organization
- cellular response to corticosteroid stimulus

Downregulated genes were enriched for processes related to:

- developmental signaling
- cell guidance
- morphogenesis
- chemotaxis
- cell migration
- cell motility

Some downregulated GO terms included nervous-system-development-like annotations, such as axon development or axon guidance. These should not be interpreted literally as neuronal differentiation in airway smooth muscle cells. A more cautious interpretation is that dexamethasone downregulated genes involved in cell guidance, migration, morphogenesis, and developmental signaling.

Main figures:

- `results/figures/task12_go_bp_up_dotplot.pdf`
- `results/figures/task12_go_bp_down_dotplot.pdf`

## 12. GSEA

GSEA was performed using the full ranked gene list based on the DESeq2 Wald statistic.

Summary:

| Item | Value |
|---|---:|
| Ranked genes used for GSEA | 16571 |
| Significant GSEA GO BP terms | 229 |
| Positive GSEA terms | 185 |
| Negative GSEA terms | 44 |

GSEA supported activation of processes related to:

- hormone-response-related biological processes
- glucose metabolism
- D-glucose import
- carbohydrate metabolism
- energy metabolism
- cell-substrate adhesion
- muscle system process
- response to oxidative stress

Some negative GSEA terms involved developmental or morphogenesis-related processes.

GO BP GSEA was useful as a global trend check, but its interpretation requires caution because GO terms are often broad, redundant, or context-dependent.

Main figures:

- `results/figures/task12_gsea_go_bp_dotplot.pdf`
- `results/figures/task12_gsea_top_positive_enrichment_plot.pdf`
- `results/figures/task12_gsea_top_negative_enrichment_plot.pdf`

## 13. Biological Interpretation

Overall, dexamethasone treatment caused broad transcriptional remodeling in airway smooth muscle cells.

The strongest signals included:

1. Activation of glucocorticoid-responsive genes  
   Examples include ZBTB16, DUSP1, PER1, and CRISPLD2.

2. Changes in extracellular matrix and adhesion-related processes  
   GO enrichment highlighted cell-substrate adhesion and extracellular matrix organization.

3. Changes in oxidative stress and metabolism  
   GSEA and GO enrichment suggested changes in oxidative stress response, glucose transport, carbohydrate metabolism, and energy metabolism.

4. Suppression of migration, chemotaxis, and morphogenesis-related programs  
   Downregulated genes were enriched for cell migration, guidance, morphogenesis, and developmental signaling-related terms.

5. Downregulation of VCAM1  
   VCAM1 showed strong downregulation and may reflect reduced adhesion or inflammatory signaling after dexamethasone treatment.

## 14. Candidate Genes for Experimental Validation

Potential validation candidates include:

| Gene | Suggested validation |
|---|---|
| ZBTB16 | qPCR for glucocorticoid-responsive upregulation |
| DUSP1 | qPCR for stress/MAPK-related glucocorticoid response |
| PER1 | qPCR for hormone/circadian-related response |
| CRISPLD2 | qPCR; known airway glucocorticoid-responsive candidate |
| VCAM1 | qPCR, flow cytometry, or immunofluorescence for downregulation |
| ADAMTS1 | qPCR for extracellular matrix remodeling |
| GPX3 | qPCR for oxidative stress / antioxidant response |
| SERPINA3 | qPCR for secreted or inflammatory-response candidate |

## 15. Limitations

This project has several limitations:

1. The analysis used a public dataset and did not include new experimental validation.
2. The sample size is small, with four paired cell lines.
3. GO enrichment results can be broad and redundant.
4. Some GO terms should not be interpreted literally without biological context.
5. RNA-seq measures mRNA abundance, not protein abundance.
6. Candidate genes require experimental validation before mechanistic conclusions can be made.

## 16. Future Directions

Possible next steps include:

1. Validate selected genes by qPCR.
2. Check whether protein-level changes match RNA-level changes.
3. Perform Hallmark or Reactome GSEA for cleaner pathway interpretation.
4. Compare results with the original airway dexamethasone study.
5. Analyze related public datasets to test whether the candidate genes are reproducible.
6. Extend the workflow to single-cell RNA-seq or disease-relevant airway datasets.

## 17. Project Summary

This project completed a full beginner-level bulk RNA-seq analysis workflow from raw FASTQ files to biological interpretation.

The analysis showed that dexamethasone treatment strongly affects gene expression in airway smooth muscle cells, inducing glucocorticoid-responsive genes and altering extracellular matrix, adhesion, oxidative stress, metabolism, and migration-related biological programs.

The final outputs include quality control reports, quantification summaries, gene-level expression matrices, DESeq2 results, visualization figures, enrichment results, and candidate genes for potential validation.
