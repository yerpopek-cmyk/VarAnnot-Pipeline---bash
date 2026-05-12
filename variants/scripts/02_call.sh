#!/usr/bin/env bash
# 02_call.sh — Variant calling with FreeBayes
# Output: ${RUN_DIR}/2_vcf_raw/${SAMPLE_ID}.raw.vcf.gz
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../workflow/config.sh"

BAM_IN="${BAM_INPUT:-${RUN_DIR}/1_bams/${SAMPLE_ID}.bam}"
VCF_OUT="${RUN_DIR}/2_vcf_raw/${SAMPLE_ID}.raw.vcf.gz"

echo "[02] FreeBayes variant calling..."
freebayes \
    -f "$REF_FASTA" \
    -p 2 \
    --min-base-quality 20 \
    --min-alternate-count 2 \
    --min-alternate-fraction "$MIN_AB" \
    "$BAM_IN" \
  | bcftools sort --max-mem 1G \
  | bgzip -@ "$THREADS" > "$VCF_OUT"

tabix -p vcf "$VCF_OUT"
echo "[02] Variants found: $(bcftools view -H "$VCF_OUT" | wc -l)"
echo "[02] Raw VCF: $VCF_OUT"
