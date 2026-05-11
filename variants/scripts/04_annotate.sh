#!/usr/bin/env bash
# 04_annotate.sh — функциональная аннотация вариантов (Ensembl VEP)
#
# VEP в одном запуске заменяет:
#   snpEff + bcftools annotate gnomAD + bcftools annotate ClinVar
#
# Выход: ${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz
#
# Примечание: нужен локальный кэш VEP (~15 ГБ для GRCh38).
# Установка кэша:
#   vep_install -a cf -s homo_sapiens -y GRCh38 \
#               --CACHEDIR data/db/vep_cache --NO_HTSLIB
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../workflow/config.sh"

VCF_IN="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.filtered.vcf.gz"
VCF_OUT="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz"

# --- Проверка кэша ---
if [[ ! -d "$VEP_CACHE_DIR" ]]; then
    echo "⚠️ Кэш VEP не найден: $VEP_CACHE_DIR"
    echo "   Начинаю автоматическую загрузку кэша (это может занять время)..."
    mkdir -p "$VEP_CACHE_DIR"
    vep_install -a cf -s homo_sapiens -y "${VEP_ASSEMBLY}" \
                --CACHEDIR "${VEP_CACHE_DIR}" --NO_HTSLIB
    echo "✅ Кэш успешно установлен!"
fi

echo "[04] VEP: аннотация вариантов (последствия + gnomAD AF + ClinVar)..."
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
echo "[04] Аннотированный VCF: $VCF_OUT"
