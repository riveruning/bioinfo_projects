# RNA-seq Foundation Project: Airway Dexamethasone Response

This repository contains a beginner-friendly but complete bulk RNA-seq analysis project.

The project analyzes public RNA-seq data from human airway smooth muscle cells treated with dexamethasone and reproduces a standard RNA-seq workflow from raw FASTQ files to biological interpretation.

## Project Goal

The main question is:

> Which genes and biological processes are affected by dexamethasone treatment in airway smooth muscle cells?

## Workflow

~~~text
FASTQ
↓
FastQC / MultiQC
↓
Salmon
↓
tximport
↓
DESeq2
↓
Visualization
↓
GO enrichment / GSEA
↓
Biological interpretation
~~~

## Dataset

Eight paired-end RNA-seq samples were analyzed.

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

The sample metadata is stored in:

- `metadata/samples.csv`

## Repository Structure

~~~text
metadata/         sample metadata and gene annotation files
scripts/          analysis scripts
logs/             command-line logs
notes/            task-by-task learning notes
results/qc/       FastQC and MultiQC results
results/tximport/ gene-level expression matrices
results/deseq2/   DESeq2 differential expression results
results/task11/   candidate gene tables and visualization summaries
results/task12/   GO enrichment and GSEA results
results/task13/   biological interpretation summaries
results/figures/  PDF figures
reports/          final project report
~~~

Large raw data and reference files are stored outside the repository and linked locally.

## Main Scripts

| Script | Purpose |
|---|---|
| `scripts/download_airway_fastq.sh` | download FASTQ files |
| `scripts/download_failed_fastq_ena.sh` | retry failed downloads using ENA |
| `scripts/run_fastqc_multiqc_all.sh` | run FastQC and MultiQC |
| `scripts/run_salmon_quant_all.sh` | run Salmon quantification |
| `scripts/run_tximport_gene_matrix.R` | import Salmon results with tximport |
| `scripts/run_deseq2_dex_vs_untrt.R` | run DESeq2 differential expression analysis |
| `scripts/visualize_deseq2_results_task11.R` | generate volcano plot, heatmap, and selected gene plots |
| `scripts/run_go_gsea_task12.R` | run GO enrichment and GSEA |
| `scripts/summarize_biology_task13.R` | summarize key genes, pathways, and validation candidates |

## Key Results

### Salmon Quantification

All eight samples showed high Salmon mapping rates, approximately 94.7% to 96.5%.

The mapping rate summary is stored in:

- `results/salmon_summary/all_samples_mapping_rate.tsv`

### tximport

tximport summarized transcript-level Salmon results to gene-level matrices.

Summary:

- samples: 8
- genes in count matrix: 82323
- transcripts in quant.sf: 654828

Main files:

- `results/tximport/gene_counts.tsv`
- `results/tximport/gene_tpm.tsv`
- `results/tximport/txi_salmon_gene_level.rds`

### DESeq2

DESeq2 was used with the design:

~~~r
design = ~ cell_line + dex
~~~

This controls for baseline differences between cell lines and estimates the dexamethasone treatment effect.

Summary:

| Item | Value |
|---|---:|
| Genes after filtering | 20407 |
| Significant genes, padj < 0.05 | 3711 |
| Upregulated significant genes | 2077 |
| Downregulated significant genes | 1634 |
| Significant genes with padj < 0.05 and abs(log2FC) >= 1 | 988 |

Main result file:

- `results/deseq2/deseq2_results_dex_trt_vs_untrt.tsv`

## Key Figures

Main figures include:

- `results/figures/task10_pca_plot.pdf`
- `results/figures/task10_ma_plot.pdf`
- `results/figures/task10_sample_distance_heatmap.pdf`
- `results/figures/task11_volcano_plot.pdf`
- `results/figures/task11_top50_candidate_gene_heatmap.pdf`
- `results/figures/task11_selected_gene_expression.pdf`
- `results/figures/task12_go_bp_up_dotplot.pdf`
- `results/figures/task12_go_bp_down_dotplot.pdf`
- `results/figures/task12_gsea_go_bp_dotplot.pdf`

## Key Biological Findings

Dexamethasone treatment induced strong transcriptional changes in airway smooth muscle cells.

Important upregulated genes included:

- ZBTB16
- DUSP1
- PER1
- CRISPLD2
- ADAMTS1
- GPX3
- SERPINA3

Important downregulated genes included:

- VCAM1

GO enrichment and GSEA suggested that dexamethasone treatment affected:

- glucocorticoid/corticosteroid response
- extracellular matrix organization
- cell-substrate adhesion
- muscle-related biological processes
- oxidative stress response
- glucose and energy metabolism
- migration, chemotaxis, and morphogenesis-related programs

## Candidate Genes for Validation

Potential experimental validation candidates include:

| Gene | Suggested validation |
|---|---|
| ZBTB16 | qPCR for glucocorticoid-responsive upregulation |
| DUSP1 | qPCR for stress/MAPK-related response |
| PER1 | qPCR for hormone/circadian response |
| CRISPLD2 | qPCR; known airway glucocorticoid-responsive candidate |
| VCAM1 | qPCR, flow cytometry, or immunofluorescence |
| ADAMTS1 | qPCR for extracellular matrix remodeling |
| GPX3 | qPCR for oxidative stress response |
| SERPINA3 | qPCR for secreted/inflammatory-response candidate |

## Final Report

The full project report is available at:

- `reports/airway_dexamethasone_rnaseq_report.md`

## Notes

Task-by-task learning notes are stored in:

- `notes/`

These notes document the reasoning, commands, issues, fixes, and interpretation for each step of the project.

## Limitations

This project is an analysis of a public dataset and does not include new experimental validation.

GO enrichment terms should be interpreted carefully because some terms are broad, redundant, or context-dependent.

RNA-seq measures mRNA abundance, not protein abundance.

## Future Work

Possible next steps include:

- qPCR validation of selected candidate genes
- protein-level validation
- Hallmark or Reactome GSEA
- comparison with the original airway dexamethasone study
- extension to single-cell or disease-relevant airway datasets
