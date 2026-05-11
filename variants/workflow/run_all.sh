#!/usr/bin/env bash
# =============================================================================
# variants/workflow/run_all.sh — оркестратор пайплайна
#
# Использование:
#   bash run_all.sh [--from STEP] [--help]
#
# Шаги: 01 02 03 04 05
# Пример перезапуска с шага 03:
#   RUN_DIR=outputs/run_20240501_120000 bash run_all.sh --from 03
# =============================================================================
set -euo pipefail

WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VARIANTS_DIR="$(dirname "$WORKFLOW_DIR")"
CONFIG="${WORKFLOW_DIR}/config.sh"

# --- Проверка конфига ---
if [[ ! -f "$CONFIG" ]]; then
    echo "❌ config.sh не найден. Скопируй шаблон:"
    echo "   cp ${WORKFLOW_DIR}/config.sh.example ${WORKFLOW_DIR}/config.sh"
    exit 1
fi
source "$CONFIG"

# --- Аргументы ---
FROM_STEP="01"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM_STEP="$2"; shift 2 ;;
        --from=*) FROM_STEP="${1#*=}"; shift ;;
        -h|--help)
            sed -n '/^# Использование/,/^# ===/p' "$0" | head -n 8; exit 0 ;;
        *) echo "❌ Неизвестный аргумент: $1"; exit 1 ;;
    esac
done

# --- Проверка входных данных ---
if [[ "$FROM_STEP" == "01" ]]; then
    if [[ -n "${BAM_INPUT:-}" ]]; then
        if [[ ! -f "$BAM_INPUT" ]]; then
            echo "❌ Ошибка: Указан BAM_INPUT, но файл не найден: $BAM_INPUT"
            exit 1
        fi
    else
        if [[ ! -f "$READS_R1" || ! -f "$READS_R2" ]]; then
            echo "❌ Ошибка: Входные FASTQ файлы не найдены:"
            echo "   R1: $READS_R1"
            echo "   R2: $READS_R2"
            exit 1
        fi
    fi
fi

if [[ ! -f "$REF_FASTA" ]]; then
    echo "❌ Ошибка: Референсный геном не найден: $REF_FASTA"
    exit 1
fi

# --- RUN_DIR: новый или переданный извне ---
if [[ -z "${RUN_DIR:-}" ]]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    export RUN_DIR="${VARIANTS_DIR}/outputs/run_${TIMESTAMP}"
fi
mkdir -p "${RUN_DIR}"/{1_bams,2_vcf_raw,3_vcf_filtered,4_reports}

echo "============================================"
echo " Variant Pipeline"
echo " Образец : $SAMPLE_ID"
echo " Папка   : $RUN_DIR"
echo " Потоков : $THREADS"
echo " С шага  : $FROM_STEP"
echo "============================================"

_run() {
    local step="$1"; local label="$2"
    if [[ "$step" < "$FROM_STEP" ]]; then
        echo "[SKIP] Шаг $step ($label)"
        return
    fi
    echo ""
    echo ">>> Шаг $step — $label"
    bash "${VARIANTS_DIR}/scripts/${step}_${label}.sh" \
        2>&1 | tee "${RUN_DIR}/step${step}.log"
}

_run 01 align
_run 02 call
_run 03 filter
_run 04 annotate

# Шаг 05 — Python
if [[ ! "05" < "$FROM_STEP" ]]; then
    echo ""
    echo ">>> Шаг 05 — report"
    python3 "${VARIANTS_DIR}/scripts/05_report.py" \
        "${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz" \
        "${RUN_DIR}/4_reports" \
        --top "${TOP_VARIANTS:-20}" \
        2>&1 | tee "${RUN_DIR}/step05.log"
fi

echo ""
echo "============================================"
echo "✅ Готово! Результаты:"
ls -lh "${RUN_DIR}/4_reports/"
echo "============================================"
