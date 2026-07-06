# Task 02: Real FASTQ Data

## 1. Task Goal

The goal of this task is to obtain, store, and inspect real RNA-seq FASTQ data from a public dataset.

## 2. Dataset

- Dataset: GSE52778
- Run accession: SRR1039508
- Organism: Homo sapiens
- Data type: bulk RNA-seq
- Sequencing layout: paired-end

## 3. Key Terms

- FASTQ: 保存测序读段和质量值的文件格式。
- read: 测序仪读出来的一小段序列，中文可理解为“测序读段”。
- FASTQ.gz: gzip 压缩后的 FASTQ 文件。
- paired-end: 双端测序，从同一个片段两端分别测序。
- R1: paired-end 中的第一端 reads。
- R2: paired-end 中的第二端 reads。
- metadata: 样本信息表，用来记录样本编号、分组、处理条件等。
- symbolic link: 软链接，类似快捷方式，可以让项目目录指向 D 盘中的真实数据目录。

## 4. Data Storage Strategy

Large FASTQ files are stored on the D drive instead of inside the WSL project directory.

The symbolic link is:

rnaseq_foundation/data/raw_fastq -> /mnt/d/bioinfo_data/rnaseq_foundation/raw_fastq

This means I can access FASTQ files from the project path:

data/raw_fastq/

but the actual files are stored on the D drive.

## 5. Downloaded Files

The complete paired-end FASTQ.gz files were downloaded:

- SRR1039508_1.fastq.gz
- SRR1039508_2.fastq.gz

File sizes:

- SRR1039508_1.fastq.gz: about 1.2G
- SRR1039508_2.fastq.gz: about 1.2G
- Total: about 2.4G

## 6. Inspecting FASTQ Content

Command used:

gzip -dc data/raw_fastq/SRR1039508_1.fastq.gz | head -n 8

This command temporarily decompresses the FASTQ.gz file and prints the first 8 lines.

Because one FASTQ read contains 4 lines, 8 lines correspond to 2 reads.

## 7. FASTQ Four-line Structure

Each read contains four lines:

1. read ID: the identifier of the sequencing read
2. sequence: the nucleotide sequence
3. separator: usually "+"
4. quality score: encoded base quality values

Example structure:

@read_id
ACTG...
+
IIII...

## 8. Important Observation

Some reads may contain N.

N means the sequencing machine could not confidently determine whether the base was A, T, C, or G.

Quality score characters such as #, H, I, and J encode base quality values.

In the observed FASTQ output, positions with N often correspond to low quality score characters such as #.

## 9. File Integrity Check

The downloaded FASTQ.gz files were checked using:

gzip -t data/raw_fastq/SRR1039508_1.fastq.gz && echo "R1 OK"
gzip -t data/raw_fastq/SRR1039508_2.fastq.gz && echo "R2 OK"

Both files passed the gzip integrity check.

## 10. Core Conclusion

This task successfully replaced toy FASTQ files with real public RNA-seq data.

The project now has a complete paired-end FASTQ input for one real sequencing run.

The next step is to run FastQC to evaluate the quality of these raw sequencing files.
