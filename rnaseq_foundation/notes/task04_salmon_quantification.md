# Task 04: Salmon Transcript Quantification

## 1. Task Goal

The goal of this task is to quantify transcript-level expression from paired-end RNA-seq FASTQ files using Salmon.

Input files:

- data/raw_fastq/SRR1039508_1.fastq.gz
- data/raw_fastq/SRR1039508_2.fastq.gz

Reference:

- GENCODE human transcript FASTA
- Salmon transcriptome index

Main output:

- results/salmon/SRR1039508/quant.sf

Summary outputs:

- results/salmon_summary/SRR1039508_top20_by_TPM.tsv
- results/salmon_summary/SRR1039508_mapping_rate.txt

## 2. Key Terms

- Salmon: a lightweight RNA-seq transcript quantification tool.
- transcript: a transcribed RNA isoform from a gene.
- transcriptome: the collection of transcript sequences.
- index: a pre-built searchable structure used to speed up read mapping.
- k-mer: a short sequence of length k.
- quantification: estimating expression abundance.
- TPM: Transcripts Per Million, a normalized expression unit.
- NumReads: estimated number of reads assigned to a transcript.
- mapping rate: the percentage of reads that can be matched to the reference transcriptome.
- gene-level expression: expression summarized at the gene level.
- transcript-level expression: expression estimated for each transcript isoform.

## 3. Build Salmon Index

Command used:

salmon index \
  -t reference/gencode.v50.transcripts.fa.gz \
  -i salmon_index/gencode_v50_salmon_index \
  -k 31

The Salmon index was built from the GENCODE human transcriptome reference.

The index was successfully built with 654828 references.

## 4. Run Salmon Quantification

Command used:

salmon quant \
  -i salmon_index/gencode_v50_salmon_index \
  -l A \
  -1 data/raw_fastq/SRR1039508_1.fastq.gz \
  -2 data/raw_fastq/SRR1039508_2.fastq.gz \
  -p 4 \
  --validateMappings \
  -o results/salmon/SRR1039508

## 5. Main Result

The Salmon quantification generated:

- quant.sf
- cmd_info.json
- lib_format_counts.json
- logs/
- aux_info/

The most important output file is quant.sf.

The columns in quant.sf are:

- Name: transcript ID and annotation information
- Length: original transcript length
- EffectiveLength: length corrected for fragment-length effects
- TPM: normalized transcript expression
- NumReads: estimated reads assigned to the transcript

## 6. Mapping Rate

The mapping rate was:

mapping rate: 95.4236%

This indicates that most reads were successfully matched to the human reference transcriptome.

## 7. Top Expressed Transcripts

The top transcripts by TPM included:

- FTL
- MT-TC
- EEF1A1
- MT-TY
- FTH1
- MT-CO1
- MT-TN
- MT-ATP8
- S100A6
- MT-CO3
- MT-ATP6

Many of these are common high-expression genes or mitochondrial transcripts in RNA-seq data.

## 8. Interpretation

The high mapping rate suggests that the FASTQ files, paired-end read pairing, and human transcriptome reference are compatible.

The quant.sf file provides transcript-level expression estimates.

For downstream differential expression analysis, transcript-level Salmon results can later be imported and summarized to gene-level counts using tximport and DESeq2.

## 9. Conclusion

Salmon quantification was completed successfully for SRR1039508.

This task converted raw FASTQ reads into transcript-level expression estimates.
