# Task 08: Batch Salmon Quantification for All Airway Samples

## 1. Task Goal

The goal of this task was to run Salmon quantification for all eight airway RNA-seq samples.

Each sample has paired-end FASTQ files:

- sample_1.fastq.gz
- sample_2.fastq.gz

Salmon was used to estimate transcript-level abundance for each sample.

## 2. Input Files

Input FASTQ files were stored under:

- data/raw_fastq/

This is a symbolic link to:

- /mnt/d/bioinfo_data/rnaseq_foundation/raw_fastq

The Salmon index used was:

- salmon_index/gencode_v50_salmon_index

The sample list was read from:

- metadata/srr_ids.txt

## 3. Output Files

Salmon output directories were generated under:

- results/salmon/

Each sample has a quantification result:

- results/salmon/SRR*/quant.sf

A mapping rate summary table was generated:

- results/salmon_summary/all_samples_mapping_rate.tsv

## 4. Mapping Rate Summary

All eight samples showed high mapping rates:

| sample_id | percent_mapped |
|---|---:|
| SRR1039508 | 95.42 |
| SRR1039509 | 94.70 |
| SRR1039512 | 96.10 |
| SRR1039513 | 96.54 |
| SRR1039516 | 95.62 |
| SRR1039517 | 96.02 |
| SRR1039520 | 96.03 |
| SRR1039521 | 96.28 |

These results suggest that the FASTQ files, paired-end structure, and transcriptome reference are compatible.

## 5. Issue Encountered

During the first batch run, SRR1039513 failed with a paired-end record size error.

The original SRR1039513 FASTQ files were likely not properly paired, even though gzip integrity checks passed.

To fix this, SRR1039513 was re-downloaded from ENA and Salmon was rerun.

The rerun completed successfully with a mapping rate of 96.54%.

## 6. Important Note

The full Salmon output directories are large and should not be committed to GitHub.

The GitHub repository should only track:

- scripts
- logs
- summary tables
- notes

## 7. Next Step

The next step is Task 09: use tximport to import Salmon quant.sf files and generate a gene-level count matrix for DESeq2.
