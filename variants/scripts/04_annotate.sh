#!/usr/bin/env bash
# 04_annotate.sh — functional variant annotation (Ensembl VEP)
#
# VEP replaces the following in one run:
#   snpEff + bcftools annotate gnomAD + bcftools annotate ClinVar
#
# Output: ${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz
#
# Note: requires local VEP cache (~15 GB for GRCh38).
# Cache installation:
#   vep_install -a cf -s homo_sapiens -y GRCh38 \
#               --CACHEDIR data/db/vep_cache --NO_HTSLIB
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../workflow/config.sh"

VCF_IN="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.filtered.vcf.gz"
VCF_OUT="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz"

# --- Cache Check ---
if [[ ! -d "$VEP_CACHE_DIR" ]]; then
    echo "⚠️ VEP cache not found: $VEP_CACHE_DIR"
    echo "   Starting automatic cache download (this may take a while)..."
    mkdir -p "$VEP_CACHE_DIR"
    vep_install -a cf -s homo_sapiens -y "${VEP_ASSEMBLY}" \
                --CACHEDIR "${VEP_CACHE_DIR}" --NO_HTSLIB
    echo "✅ Cache successfully installed!"
fi

echo "[04] VEP: variant annotation (consequences + gnomAD AF + ClinVar)..."
vep \
    --input_file  "$VCF_IN" \
    --output_file "$VCF_OUT" \
    --format vcf \
    --vcf \
    --compress_output bgzip \
    --offline \
    --cache \
    --dir_cache "$VEP_CACHE_DIR" \
    --assembly  "$VEP_ASSEMBLY" \
    --fork      "$THREADS" \
    --canonical \
    --hgvs \
    --check_existing \
    --af_gnomade \
    --af_gnomadg \
    --no_stats \
    --quiet \
    --force_overwrite \
    --synonyms "data/db/chr_synonyms.txt" \
    --fields "Consequence,IMPACT,SYMBOL,Gene,Feature,HGVSc,HGVSp,gnomADe_AF,gnomADg_AF,CLIN_SIG"

tabix -p vcf "$VCF_OUT"
echo "[04] Annotated VCF: $VCF_OUT"
