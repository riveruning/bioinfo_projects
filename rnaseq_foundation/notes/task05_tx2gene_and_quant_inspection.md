# Task 05: tx2gene Table and Salmon Quant Inspection

## 1. Task Goal

The goal of this task is to prepare the bridge between Salmon transcript-level quantification and downstream gene-level analysis.

Salmon outputs transcript-level expression estimates in quant.sf.

Downstream differential expression analysis tools such as DESeq2 usually work with gene-level count matrices.

Therefore, a tx2gene table is needed to map transcript IDs to gene IDs.

## 2. Key Terms

- transcript: a transcript isoform produced from a gene.
- gene: a genomic unit that may produce one or more transcripts.
- tx2gene: a table mapping transcript IDs to gene IDs.
- tximport: an R package used to import transcript-level quantification and summarize it to the gene level.
- transcript-level expression: expression estimated for each transcript isoform.
- gene-level expression: expression summarized for each gene.
- TPM: Transcripts Per Million, a normalized expression value.
- NumReads: Salmon's estimated number of reads assigned to a transcript.
- Rscript: a command-line tool for running R scripts.

## 3. Files Generated

The following files were generated:

- metadata/tx2gene_full.tsv
- metadata/tx2gene.tsv
- scripts/inspect_salmon_quant.R

metadata/tx2gene_full.tsv contains:

- transcript_id
- gene_id
- gene_name
- gene_type

metadata/tx2gene.tsv contains:

- transcript_id
- gene_id

## 4. R Quantification Inspection

The Salmon quantification result was inspected using:

Rscript scripts/inspect_salmon_quant.R

The output showed:

- Number of transcripts: 654828
- Total estimated NumReads: 21885909
- Detected transcripts with TPM > 0: 128858

The top transcripts by TPM included FTL, MT-TC, EEF1A1, MT-TY, FTH1, MT-CO1, MT-ATP8, S100A6, MT-CO3, MT-ATP6, RPS27, RPL41, RPS18, TMSB4X, and VIM.

## 5. Interpretation

The tx2gene table links transcript-level Salmon results to gene IDs.

This mapping will be required by tximport to summarize transcript-level quantification into gene-level expression.

The R inspection script confirms that quant.sf can be successfully loaded and analyzed in R.

## 6. Conclusion

This task prepared the project for downstream gene-level RNA-seq analysis.

The next major step is to quantify multiple samples and then use tximport and DESeq2 for differential expression analysis.
