# Reproducibility Notes

This file summarizes how the project can be reproduced from the repository structure.

## 1. Project Root

Run commands from:

~~~bash
cd ~/bioinfo_projects/rnaseq_foundation
~~~

## 2. Main Workflow

The analysis workflow is organized into scripts:

~~~text
1. Download FASTQ files
2. Run FastQC / MultiQC
3. Run Salmon quantification
4. Import Salmon results with tximport
5. Run DESeq2 differential expression analysis
6. Visualize DESeq2 results
7. Run GO enrichment and GSEA
8. Summarize biological interpretation
~~~

## 3. Main Scripts

~~~bash
scripts/download_airway_fastq.sh
scripts/download_failed_fastq_ena.sh
scripts/run_fastqc_multiqc_all.sh
scripts/run_salmon_quant_all.sh
scripts/run_tximport_gene_matrix.R
scripts/run_deseq2_dex_vs_untrt.R
scripts/visualize_deseq2_results_task11.R
scripts/run_go_gsea_task12.R
scripts/summarize_biology_task13.R
~~~

## 4. Main Inputs

Small metadata and annotation files are stored in:

~~~text
metadata/samples.csv
metadata/srr_ids.txt
metadata/tx2gene.tsv
metadata/tx2gene_full.tsv
~~~

Large raw data and reference files are not stored in Git.

They are expected to be available locally through these paths or symbolic links:

~~~text
data/raw_fastq/
reference/
salmon_index/
~~~

## 5. Main Outputs

Important results are stored in:

~~~text
results/qc/
results/salmon_summary/
results/tximport/
results/deseq2/
results/task11/
results/task12/
results/task13/
results/figures/
reports/
~~~

## 6. Important Notes

- Raw FASTQ files are not tracked by Git.
- GENCODE reference files are not tracked by Git.
- Salmon index files are not tracked by Git.
- Full Salmon output directories are not tracked by Git.
- Summary tables, figures, logs, notes, and reports are tracked.

## 7. Final Report

The final project report is:

~~~text
reports/airway_dexamethasone_rnaseq_report.md
~~~
