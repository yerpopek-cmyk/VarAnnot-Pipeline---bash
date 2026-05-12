#!/usr/bin/env bash
# 01_align.sh — BWA: alignment → MarkDup → index
# Output BAM: ${RUN_DIR}/1_bams/${SAMPLE_ID}.bam
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../workflow/config.sh"

OUT_BAM="${RUN_DIR}/1_bams/${SAMPLE_ID}.bam"

# --- Skip step if a ready BAM is used ---
if [[ -n "${BAM_INPUT:-}" ]]; then
    echo "[01] BAM_INPUT is set ($BAM_INPUT). Skipping FASTQ alignment."
    exit 0
fi

# --- Indexing the reference (if necessary) ---
[[ -f "${REF_FASTA}.fai" ]]              || samtools faidx "$REF_FASTA"
[[ -f "${REF_FASTA}.bwt" ]]              || bwa index "$REF_FASTA"

echo "[01] BWA (classic) → sort → markdup..."
bwa mem \
    -t "$THREADS" \
    -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1" \
    "$REF_FASTA" "$READS_R1" "$READS_R2" \
  | samtools fixmate -m -u - - \
  | samtools sort -@ "$THREADS" -m "$MEM_SORT" \
  | samtools markdup -@ "$THREADS" --write-index - "$OUT_BAM"

echo "[01] Alignment flags:"
samtools flagstat "$OUT_BAM" | tee "${RUN_DIR}/1_bams/${SAMPLE_ID}.flagstat"
echo "[01] BAM is ready: $OUT_BAM"
