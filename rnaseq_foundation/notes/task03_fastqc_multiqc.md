# Task 03: FastQC and MultiQC

## 1. Task Goal

The goal of this task is to perform quality control on real paired-end RNA-seq FASTQ files.

Input files:

- data/raw_fastq/SRR1039508_1.fastq.gz
- data/raw_fastq/SRR1039508_2.fastq.gz

Output files:

- results/qc/fastqc/SRR1039508_1_fastqc.html
- results/qc/fastqc/SRR1039508_1_fastqc.zip
- results/qc/fastqc/SRR1039508_2_fastqc.html
- results/qc/fastqc/SRR1039508_2_fastqc.zip
- results/qc/multiqc/multiqc_report.html

## 2. Key Terms

- FastQC: FASTQ quality control tool.
- MultiQC: tool for summarizing multiple bioinformatics reports into one report.
- QC: quality control.
- read: a sequencing read.
- R1: the first read file in paired-end sequencing.
- R2: the second read file in paired-end sequencing.
- adapter: artificial sequence added during library preparation.
- duplication: repeated sequencing reads.
- overrepresented sequence: a sequence that appears more frequently than expected.
- GC content: the percentage of G and C bases in reads.
- N content: the percentage of unknown bases.

## 3. Commands Used

Create output directories:

mkdir -p results/qc/fastqc
mkdir -p results/qc/multiqc

Run FastQC:

fastqc data/raw_fastq/SRR1039508_1.fastq.gz \
       data/raw_fastq/SRR1039508_2.fastq.gz \
       -o results/qc/fastqc

Run MultiQC:

multiqc results/qc/fastqc -o results/qc/multiqc

## 4. Output Files

FastQC generated one HTML report and one ZIP result file for each FASTQ file.

MultiQC generated:

- multiqc_report.html
- multiqc_data/

## 5. Main QC Findings

The paired-end FASTQ files were successfully processed by FastQC and MultiQC.

Main observations:

1. The overall sequencing quality is acceptable.
2. Per-base sequence quality passed for both R1 and R2.
3. Per-sequence quality scores passed for both R1 and R2.
4. Per-sequence GC content passed for both R1 and R2.
5. Adapter content passed, suggesting no strong adapter contamination.
6. Sequence duplication levels were elevated.
7. Overrepresented sequences were detected, but the total proportion was low.
8. R2 showed slightly more QC warnings than R1, which is common in paired-end sequencing.

## 6. Interpretation

The high duplication level should be recorded, but it does not necessarily mean the RNA-seq data is unusable.

In RNA-seq, highly expressed transcripts can naturally produce many repeated reads. Therefore, duplication warnings should be interpreted together with the biological context.

Adapter content passed, so there is no strong need to perform adapter trimming at this stage.

## 7. Conclusion

The FASTQ files passed the main quality checks and are suitable for downstream RNA-seq analysis.

The next step is expression quantification using Salmon.
