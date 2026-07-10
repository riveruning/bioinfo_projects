# Task 07: FastQC and MultiQC for All Airway Samples

## 1. Task Goal

The goal of this task was to run quality control for all eight airway RNA-seq samples.

Each sample has paired-end FASTQ files:

- R1: sample_1.fastq.gz
- R2: sample_2.fastq.gz

Therefore, 16 FASTQ files were checked in total.

## 2. Tools Used

- FastQC: quality control for individual FASTQ files
- MultiQC: summary report for all FastQC results

## 3. Input Files

Input FASTQ files were stored outside the Git repository under:

- /mnt/d/bioinfo_data/rnaseq_foundation/raw_fastq

The project uses a symbolic link:

- data/raw_fastq

## 4. Output Files

FastQC output:

- results/qc/fastqc_all_samples/

MultiQC output:

- results/qc/multiqc_all_samples/multiqc_report.html
- results/qc/multiqc_all_samples/multiqc_report_data/

## 5. Run Summary

FastQC completed for all 16 paired-end FASTQ files.

MultiQC found 16 FastQC reports and generated a project-level QC report.

## 6. Initial Interpretation

The MultiQC report should be inspected for:

- per base sequence quality
- per sequence quality scores
- per base N content
- adapter content
- sequence duplication levels
- overrepresented sequences
- abnormal samples or outliers

If there are no severe quality problems or obvious outlier samples, the project can proceed to batch Salmon quantification.

## 7. Next Step

The next step is Task 08: run Salmon quantification for all eight samples.
