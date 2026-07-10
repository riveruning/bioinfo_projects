#!/usr/bin/env bash
set -euo pipefail

RAW_DIR=$(readlink -f data/raw_fastq)

FASTQC_OUT="results/qc/fastqc_all_samples"
MULTIQC_OUT="results/qc/multiqc_all_samples"

mkdir -p "$FASTQC_OUT"
mkdir -p "$MULTIQC_OUT"
mkdir -p logs

echo "FASTQ directory:"
echo "$RAW_DIR"
echo

echo "Checking required tools..."
command -v fastqc
command -v multiqc
echo

echo "Checking FASTQ files from metadata/srr_ids.txt..."
for srr in $(cat metadata/srr_ids.txt)
do
  r1="$RAW_DIR/${srr}_1.fastq.gz"
  r2="$RAW_DIR/${srr}_2.fastq.gz"

  if [[ ! -s "$r1" ]]; then
    echo "Missing R1: $r1"
    exit 1
  fi

  if [[ ! -s "$r2" ]]; then
    echo "Missing R2: $r2"
    exit 1
  fi

  echo "$srr OK"
done

echo
echo "Running FastQC for all FASTQ files..."
fastqc \
  -t 4 \
  -o "$FASTQC_OUT" \
  "$RAW_DIR"/SRR10395*_1.fastq.gz \
  "$RAW_DIR"/SRR10395*_2.fastq.gz

echo
echo "Running MultiQC..."
multiqc \
  "$FASTQC_OUT" \
  -o "$MULTIQC_OUT" \
  --filename multiqc_report.html \
  --force

echo
echo "Task 07 finished."
echo "FastQC output: $FASTQC_OUT"
echo "MultiQC report: $MULTIQC_OUT/multiqc_report.html"
