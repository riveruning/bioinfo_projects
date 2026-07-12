# Task 09: Import Salmon Quantification with tximport

## 1. Task Goal

The goal of this task was to import Salmon transcript-level quantification results from all eight airway RNA-seq samples and summarize them to gene-level expression matrices.

The output gene-level count matrix will be used for DESeq2 differential expression analysis.

## 2. Input Files

Input Salmon quantification files:

- results/salmon/SRR1039508/quant.sf
- results/salmon/SRR1039509/quant.sf
- results/salmon/SRR1039512/quant.sf
- results/salmon/SRR1039513/quant.sf
- results/salmon/SRR1039516/quant.sf
- results/salmon/SRR1039517/quant.sf
- results/salmon/SRR1039520/quant.sf
- results/salmon/SRR1039521/quant.sf

Sample metadata:

- metadata/samples.csv

Transcript-to-gene annotation:

- metadata/tx2gene_full.tsv

## 3. Transcript ID Matching

The original Salmon quant.sf transcript names did not exactly match the transcript_id column in tx2gene_full.tsv.

Exact transcript ID overlap rate was 0%.

This happened because Salmon used GENCODE-style transcript names containing additional information separated by vertical bars.

After parsing the GENCODE-style names, the transcript-to-gene matching rate was 100%.

## 4. tximport Output

The following files were generated:

- results/tximport/gene_counts.tsv
- results/tximport/gene_tpm.tsv
- results/tximport/gene_effective_length.tsv
- results/tximport/gene_annotation.tsv
- results/tximport/txi_salmon_gene_level.rds
- results/tximport/tximport_summary.tsv

## 5. Summary

The tximport summary was:

- samples: 8
- genes in count matrix: 82323
- transcripts in quant.sf: 654828
- transcripts in tx2gene_for_salmon: 654828

This means that all eight samples were successfully imported, and Salmon transcript-level estimates were summarized to a gene-level expression matrix.

## 6. Important Note

The count values in gene_counts.tsv are estimated counts from Salmon and tximport, so they may contain decimal values.

This is normal for transcript-level quantification methods because reads can be probabilistically assigned to transcripts and then summarized to genes.

## 7. Next Step

The next step is Task 10: differential expression analysis with DESeq2.

The main comparison will be dexamethasone-treated samples versus untreated samples, while accounting for cell line differences.
