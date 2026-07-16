# Software Environment

This project was run in WSL Ubuntu using conda environments.

## Main Tools

Command-line tools:

- FastQC
- MultiQC
- Salmon
- SRA Toolkit / fastq-dump
- wget
- gzip
- Git

R / Bioconductor packages:

- tximport
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- enrichplot
- fgsea
- ggplot2
- pheatmap
- dplyr

## Notes

The main RNA-seq workflow was developed in the `bioinfo` conda environment.

During GO enrichment setup, additional dependency issues involving `jq` and `yq` were encountered and fixed.

For long-term reproducibility, it would be better to create a fresh environment specifically for this project in the future.
