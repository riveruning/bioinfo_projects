#!/usr/bin/env bash
set -euo pipefail

RAW_DIR=$(readlink -f data/raw_fastq)
INDEX_DIR="salmon_index/gencode_v50_salmon_index"
OUT_DIR="results/salmon"
SUMMARY_DIR="results/salmon_summary"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$SUMMARY_DIR" "$LOG_DIR"

echo "FASTQ directory:"
echo "$RAW_DIR"
echo

echo "Salmon index:"
echo "$INDEX_DIR"
echo

echo "Checking Salmon..."
command -v salmon
salmon --version
echo

echo "Checking input FASTQ files..."
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

  echo "$srr FASTQ OK"
done

echo
echo "Running Salmon quantification..."

for srr in $(cat metadata/srr_ids.txt)
do
  r1="$RAW_DIR/${srr}_1.fastq.gz"
  r2="$RAW_DIR/${srr}_2.fastq.gz"
  sample_out="${OUT_DIR}/${srr}"

  echo "=============================="
  echo "Processing $srr"
  echo "=============================="

  if [[ -s "${sample_out}/quant.sf" ]]; then
    echo "$srr already has quant.sf, skipping Salmon."
  else
    salmon quant \
      -i "$INDEX_DIR" \
      -l A \
      -1 "$r1" \
      -2 "$r2" \
      -p 4 \
      --validateMappings \
      -o "$sample_out"
  fi

  if [[ ! -s "${sample_out}/quant.sf" ]]; then
    echo "ERROR: quant.sf was not generated for $srr"
    exit 1
  fi

  echo "$srr quant.sf OK"
done

echo
echo "Creating mapping rate summary..."

python - <<'PY'
import json
from pathlib import Path

samples = [line.strip() for line in Path("metadata/srr_ids.txt").read_text().splitlines() if line.strip()]
out_path = Path("results/salmon_summary/all_samples_mapping_rate.tsv")
out_path.parent.mkdir(parents=True, exist_ok=True)

with out_path.open("w") as out:
    out.write("sample_id\tnum_processed\tnum_mapped\tpercent_mapped\n")
    for srr in samples:
        meta_path = Path("results/salmon") / srr / "aux_info" / "meta_info.json"
        if not meta_path.exists():
            raise FileNotFoundError(f"Missing meta_info.json for {srr}: {meta_path}")

        meta = json.loads(meta_path.read_text())

        num_processed = meta.get("num_processed", "NA")
        num_mapped = meta.get("num_mapped", "NA")
        percent_mapped = meta.get("percent_mapped", "NA")

        out.write(f"{srr}\t{num_processed}\t{num_mapped}\t{percent_mapped}\n")

print(f"Wrote {out_path}")
PY

echo
echo "Task 08 finished."
echo "Salmon output directory: $OUT_DIR"
echo "Mapping rate summary: $SUMMARY_DIR/all_samples_mapping_rate.tsv"
