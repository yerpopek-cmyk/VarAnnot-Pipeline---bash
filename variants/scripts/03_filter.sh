#!/usr/bin/env bash
# 03_filter.sh — нормализация VCF + двухэтапная фильтрация
#
# Этап 1: bcftools norm   — left-align, split multiallelic
# Этап 2: soft-filter     — QUAL, DP (помечаем, не удаляем — VEP их увидит)
# Этап 3: soft-filter AB  — allele balance для гетерозигот (FreeBayes: AO/RO)
#
# Итог: ${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.filtered.vcf.gz
#       (FILTER=PASS — прошли все пороги; LowQual/LowAB — помечены)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../workflow/config.sh"

VCF_IN="${RUN_DIR}/2_vcf_raw/${SAMPLE_ID}.raw.vcf.gz"
VCF_OUT="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.filtered.vcf.gz"
STATS="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.stats.txt"

echo "[03] Нормализация + фильтрация..."
bcftools norm \
    -f "$REF_FASTA" \
    -m -any \
    "$VCF_IN" \
  | bcftools filter \
    --soft-filter LowQual \
    --mode + \
    -e "QUAL < ${MIN_QUAL} || INFO/DP < ${MIN_DP} || INFO/DP > ${MAX_DP}" \
  | bcftools filter \
    --soft-filter LowAB \
    --mode + \
    -e 'GT[*]="het" && (FORMAT/AO[*:0] / (FORMAT/RO[*:0] + FORMAT/AO[*:0]) < '"${MIN_AB}"' || FORMAT/AO[*:0] / (FORMAT/RO[*:0] + FORMAT/AO[*:0]) > '"${MAX_AB}"')' \
    -Oz -o "$VCF_OUT"
tabix -p vcf "$VCF_OUT"

echo "[03] Статистика:"
bcftools stats "$VCF_OUT" | grep -E "^SN" | tee "$STATS"

TOTAL=$(bcftools view -H "$VCF_OUT" | wc -l)
PASS=$(bcftools view -f PASS -H "$VCF_OUT" | wc -l)
echo "[03] Всего: $TOTAL  |  PASS: $PASS"
echo "[03] Отфильтрованный VCF: $VCF_OUT"
