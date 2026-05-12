#!/usr/bin/env bash
# =============================================================================
# variants/workflow/run_all.sh — Pipeline Orchestrator
#
# Usage:
#   bash run_all.sh [--from STEP] [--help]
#
# Steps: 01 02 03 04 05
# Example to restart from step 03:
#   RUN_DIR=outputs/run_20240501_120000 bash run_all.sh --from 03
# =============================================================================
set -euo pipefail

WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VARIANTS_DIR="$(dirname "$WORKFLOW_DIR")"
CONFIG="${WORKFLOW_DIR}/config.sh"

# --- Config Check ---
if [[ ! -f "$CONFIG" ]]; then
    echo "❌ config.sh not found. Copy the template:"
    echo "   cp ${WORKFLOW_DIR}/config.sh.example ${WORKFLOW_DIR}/config.sh"
    exit 1
fi
source "$CONFIG"

# --- Arguments ---
FROM_STEP="01"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM_STEP="$2"; shift 2 ;;
        --from=*) FROM_STEP="${1#*=}"; shift ;;
        -h|--help)
            sed -n '/^# Usage/,/^# ===/p' "$0" | head -n 8; exit 0 ;;
        *) echo "❌ Unknown argument: $1"; exit 1 ;;
    esac
done

# --- Input Validation ---
if [[ "$FROM_STEP" == "01" ]]; then
    if [[ -n "${BAM_INPUT:-}" ]]; then
        if [[ ! -f "$BAM_INPUT" ]]; then
            echo "❌ Error: BAM_INPUT is specified, but the file was not found: $BAM_INPUT"
            exit 1
        fi
    else
        if [[ ! -f "$READS_R1" || ! -f "$READS_R2" ]]; then
            echo "❌ Error: Input FASTQ files not found:"
            echo "   R1: $READS_R1"
            echo "   R2: $READS_R2"
            exit 1
        fi
    fi
fi

if [[ ! -f "$REF_FASTA" ]]; then
    echo "❌ Error: Reference genome not found: $REF_FASTA"
    exit 1
fi

# --- RUN_DIR: create new or use external ---
if [[ -z "${RUN_DIR:-}" ]]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    export RUN_DIR="${VARIANTS_DIR}/outputs/run_${TIMESTAMP}"
fi
mkdir -p "${RUN_DIR}"/{1_bams,2_vcf_raw,3_vcf_filtered,4_reports}

echo "============================================"
echo " Variant Pipeline"
echo " Sample  : $SAMPLE_ID"
echo " Folder  : $RUN_DIR"
echo " Threads : $THREADS"
echo " Step    : $FROM_STEP"
echo "============================================"

_run() {
    local step="$1"; local label="$2"
    if [[ "$step" < "$FROM_STEP" ]]; then
        echo "[SKIP] Step $step ($label)"
        return
    fi
    echo ""
    echo ">>> Step $step — $label"
    bash "${VARIANTS_DIR}/scripts/${step}_${label}.sh" \
        2>&1 | tee "${RUN_DIR}/step${step}.log"
}

_run 01 align
_run 02 call
_run 03 filter
_run 04 annotate

# Step 05 — Python
if [[ ! "05" < "$FROM_STEP" ]]; then
    echo ""
    echo ">>> Step 05 — report"
    python3 "${VARIANTS_DIR}/scripts/05_report.py" \
        "${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz" \
        "${RUN_DIR}/4_reports" \
        --top "${TOP_VARIANTS:-20}" \
        2>&1 | tee "${RUN_DIR}/step05.log"
fi

echo ""
echo "============================================"
echo "✅ Done! Results:"
ls -lh "${RUN_DIR}/4_reports/"
echo "============================================"
