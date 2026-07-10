# Task 06: Download All Airway FASTQ Files

## 1. Task Goal

The goal of this task was to expand the project from one RNA-seq sample to all eight airway samples.

The airway dataset contains four cell lines, each with an untreated sample and a dexamethasone-treated sample.

## 2. Files Prepared

The following metadata files were prepared:

- metadata/samples.csv
- metadata/srr_ids.txt

metadata/samples.csv records:

- sample_id
- cell_line
- dex
- condition

metadata/srr_ids.txt contains the SRR run IDs used for downloading FASTQ files.

## 3. Download Strategy

The first sample was downloaded using fastq-dump.

During batch downloading, fastq-dump showed unstable network-related errors, including:

- Failed to call external services
- timeout errors
- incomplete FASTQ files

Therefore, failed samples were downloaded from ENA using wget with resume support.

ENA provides direct FASTQ.gz files, so these files can be used directly by Salmon after gzip integrity checks.

## 4. FASTQ Integrity Check

Each sample should have two paired-end FASTQ files:

- sample_1.fastq.gz
- sample_2.fastq.gz

All FASTQ files were checked using gzip -t.

A FASTQ file is only considered valid if gzip -t passes.

## 5. Current Status

All eight airway samples have been downloaded:

- SRR1039508
- SRR1039509
- SRR1039512
- SRR1039513
- SRR1039516
- SRR1039517
- SRR1039520
- SRR1039521

The next step is to run Salmon quantification for all eight samples.

## 6. Conclusion

This task completed the raw FASTQ preparation step for the full airway RNA-seq dataset.

The project is now ready for batch Salmon quantification.
