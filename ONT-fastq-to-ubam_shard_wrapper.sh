#!/usr/bin/env bash
set -euo pipefail

# shard_ont_fastq_to_ubam.sh
# Split a FASTQ into N contiguous shards (4-line records), convert each shard
# to uBAM with ONT-fastq-to-ubam.py in parallel, then merge with samtools cat.

usage() {
  cat <<EOF
Usage:
  $(basename "$0") -i INPUT_FASTQ[.gz] -o OUTPUT.ubam -s SAMPLE_NAME -n SHARDS -p /path/to/ONT-fastq-to-ubam.py [-- EXTRA_ARGS]

Required:
  -i  Input FASTQ (.fastq or .fastq.gz)
  -o  Output uBAM path
  -s  Sample name for @RG:ID and @RG:SM
  -n  Number of shards (parallel workers)
  -p  Path to ONT-fastq-to-ubam.py

Optional:
  Anything after -- is passed to ONT-fastq-to-ubam.py (e.g., --no-individual-tags)

Requires: samtools, coreutils 'split', awk, xargs; pigz recommended for .gz
EOF
}

in_fq="" ; out_bam="" ; sample="" ; shards="" ; conv_py=""
extra_args=()
while (( "$#" )); do
  case "$1" in
    -i) in_fq="$2"; shift 2 ;;
    -o) out_bam="$2"; shift 2 ;;
    -s) sample="$2"; shift 2 ;;
    -n) shards="$2"; shift 2 ;;
    -p) conv_py="$2"; shift 2 ;;
    --) shift; extra_args=("$@"); break ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$in_fq" || -z "$out_bam" || -z "$sample" || -z "$shards" || -z "$conv_py" ]] && { usage; exit 1; }
command -v samtools >/dev/null || { echo "samtools not found"; exit 1; }
[[ -x "$conv_py" || "$conv_py" =~ \.py$ ]] || { echo "Converter not found: $conv_py"; exit 1; }

tmpdir="$(mktemp -d -t ontubamshard.XXXXXX)"
cleanup(){ rm -rf "$tmpdir"; }
trap cleanup EXIT

# 1) Count total lines, compute lines per shard (multiple of 4)
if [[ "$in_fq" =~ \.gz$ ]]; then
  if command -v pigz >/dev/null; then
    total_lines=$(pigz -dc -- "$in_fq" | wc -l)
  else
    total_lines=$(gzip -dc -- "$in_fq" | wc -l)
  fi
else
  total_lines=$(wc -l < "$in_fq")
fi

if (( total_lines % 4 != 0 )); then
  echo "Error: FASTQ line count not divisible by 4 (${total_lines})." >&2
  exit 1
fi

reads=$(( total_lines / 4 ))
# ceil(reads / shards)
chunk_reads=$(( (reads + shards - 1) / shards ))
chunk_lines=$(( chunk_reads * 4 ))

echo "Total reads: $reads  | shards: $shards  | reads/shard ~= $chunk_reads"

# 2) Split into contiguous 4-line chunks
if [[ "$in_fq" =~ \.gz$ ]]; then
  if command -v pigz >/dev/null; then
    pigz -dc -- "$in_fq" | split -l "$chunk_lines" -d -a 3 - "$tmpdir/shard_"
  else
    gzip -dc -- "$in_fq" | split -l "$chunk_lines" -d -a 3 - "$tmpdir/shard_"
  fi
else
  split -l "$chunk_lines" -d -a 3 -- "$in_fq" "$tmpdir/shard_"
fi

# Rename split parts to .fastq for clarity
for f in "$tmpdir"/shard_*; do mv "$f" "${f}.fastq"; done

# 3) Convert each shard to uBAM in parallel
ls "$tmpdir"/shard_*.fastq | sort | \
  xargs -n1 -P "$shards" -I{} \
    bash -c 'python3 "'"$conv_py"'" "{}" "{}".ubam --sample-name "'"$sample"'" '"${extra_args[@]+"${extra_args[@]}"}"

# 4) Merge with samtools cat (preserves order by sorted shard index)
ls "$tmpdir"/shard_*.fastq.ubam | sort | xargs samtools cat -o "$out_bam"

echo "Done: $out_bam"