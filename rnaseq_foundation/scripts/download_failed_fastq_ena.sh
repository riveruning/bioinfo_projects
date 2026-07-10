#!/usr/bin/env bash
set -euo pipefail

RAW_DIR=$(readlink -f data/raw_fastq)
mkdir -p logs tmp

FAILED_LOG="logs/failed_srr_ids_ena.txt"
: > "$FAILED_LOG"

echo "FASTQ output directory:"
echo "$RAW_DIR"

while read -r srr
do
  echo "=============================="
  echo "Processing ${srr}"
  echo "=============================="

  r1="${RAW_DIR}/${srr}_1.fastq.gz"
  r2="${RAW_DIR}/${srr}_2.fastq.gz"

  if [[ -s "$r1" && -s "$r2" ]] && gzip -t "$r1" && gzip -t "$r2"; then
    echo "${srr}: already exists and passed gzip check. Skipping."
    continue
  fi

  echo "${srr}: removing old incomplete files..."
  rm -f "${RAW_DIR}/${srr}"*.fastq
  rm -f "${RAW_DIR}/${srr}"*.fastq.gz

  api_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"

  echo "${srr}: querying ENA..."
  curl -L -s "$api_url" > "tmp/${srr}_ena.tsv"

  urls=$(tail -n +2 "tmp/${srr}_ena.tsv" | cut -f2 | tr ';' '\n' | grep '_[12]\.fastq\.gz$' || true)

  if [[ -z "$urls" ]]; then
    echo "${srr}: no FASTQ URLs found from ENA."
    echo "$srr" >> "$FAILED_LOG"
    continue
  fi

  echo "${srr}: ENA FASTQ URLs:"
  echo "$urls"

  while read -r url
  do
    [[ -z "$url" ]] && continue

    filename=$(basename "$url")

    if [[ "$url" == http://* || "$url" == https://* || "$url" == ftp://* ]]; then
      dl_url="$url"
    else
      dl_url="https://${url}"
    fi

    echo "${srr}: downloading ${filename}"
    wget -c -P "$RAW_DIR" "$dl_url"
  done <<< "$urls"

  echo "${srr}: checking gzip integrity..."
  if [[ -s "$r1" && -s "$r2" ]] && gzip -t "$r1" && gzip -t "$r2"; then
    echo "${srr}: OK"
  else
    echo "${srr}: FAILED gzip check. Removing incomplete files."
    rm -f "${RAW_DIR}/${srr}"*.fastq
    rm -f "${RAW_DIR}/${srr}"*.fastq.gz
    echo "$srr" >> "$FAILED_LOG"
  fi

done < metadata/srr_ids_failed_retry.txt

echo "=============================="
echo "Finished ENA downloads"
echo "Failed SRR IDs, if any:"
cat "$FAILED_LOG"
