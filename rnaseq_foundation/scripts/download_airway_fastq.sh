#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw_fastq
mkdir -p logs

RAW_DIR=$(readlink -f data/raw_fastq)

echo "FASTQ output directory:"
echo "$RAW_DIR"

while read -r srr
do
  r1="${RAW_DIR}/${srr}_1.fastq.gz"
  r2="${RAW_DIR}/${srr}_2.fastq.gz"

  echo "=============================="
  echo "Processing ${srr}"
  echo "=============================="

  if [[ -s "$r1" && -s "$r2" ]]; then
    echo "${srr}: FASTQ files already exist. Skipping download."
  else
    echo "${srr}: downloading FASTQ files..."
    rm -f "${RAW_DIR}/${srr}"*.fastq "${RAW_DIR}/${srr}"*.fastq.gz
    fastq-dump --split-files --gzip --outdir "$RAW_DIR" "$srr"
  fi

  echo "${srr}: checking gzip integrity..."
  gzip -t "$r1"
  gzip -t "$r2"

  echo "${srr}: OK"
done < metadata/srr_ids.txt
