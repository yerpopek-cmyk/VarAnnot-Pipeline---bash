#!/usr/bin/env bash
# =============================================================================
# annotation/run_annotation.sh
# Функциональная аннотация бактериального генома:
#   Prokka (структурная) → DIAMOND blastp (гомология) → hmmscan (домены Pfam)
#
# Использование:
#   bash run_annotation.sh [--offline] [--max-proteins N] [--cpus N]
#
# Опции:
#   --offline          Пропустить скачивание, использовать кэш из data/db/
#   --max-proteins N   Ограничить число белков для hmmscan (удобно при тесте)
#   --cpus N           Число потоков (по умолчанию: nproc − 2, но не менее 1)
# =============================================================================
set -euo pipefail

# --- Пути ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
DB_DIR="${DATA_DIR}/db"
INPUT_DIR="${DATA_DIR}/input"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUT_DIR="${SCRIPT_DIR}/outputs/run_${TIMESTAMP}"

ASSEMBLY="${INPUT_DIR}/assembly.fasta"
PROTEIN_DB_FASTA="${DB_DIR}/swissprot_subset.fasta"
DIAMOND_DB="${DB_DIR}/swissprot_subset"   # без расширения
HMM_DB="${DB_DIR}/pfam_subset.hmm"

PROKKA_PREFIX="ANNOT"

# --- Умный подсчёт ядер ---
TOTAL_CPUS=$(nproc 2>/dev/null || echo 4)
DEFAULT_CPUS=$(( TOTAL_CPUS > 2 ? TOTAL_CPUS - 2 : TOTAL_CPUS ))
CPUS="$DEFAULT_CPUS"
OFFLINE=false
MAX_PROTEINS=0

# --- Аргументы ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --offline)           OFFLINE=true;         shift ;;
        --max-proteins)      MAX_PROTEINS="$2";    shift 2 ;;
        --max-proteins=*)    MAX_PROTEINS="${1#*=}"; shift ;;
        --cpus)              CPUS="$2";            shift 2 ;;
        --cpus=*)            CPUS="${1#*=}";       shift ;;
        -h|--help)
            sed -n '/^# Использование/,/^# ===/p' "$0" | head -n 8
            exit 0 ;;
        *)  echo "❌ Неизвестный аргумент: $1"; exit 1 ;;
    esac
done

echo "⚙️  Потоков: $CPUS / $TOTAL_CPUS  |  Выходная папка: $OUT_DIR"
mkdir -p "$DB_DIR" "$INPUT_DIR" "$OUT_DIR"/{prokka,diamond,hmmer}

# =============================================================================
# 1. Проверка инструментов
# =============================================================================
echo -e "\n🔍 Проверка инструментов..."
MISSING=()
for tool in prokka diamond hmmscan hmmpress seqkit wget bgzip; do
    command -v "$tool" &>/dev/null || MISSING+=("$tool")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "❌ Не найдены: ${MISSING[*]}"
    echo "   Активируй окружение: conda activate annotation"
    exit 1
fi
echo "✅ Все инструменты найдены"

# =============================================================================
# 2. Загрузка данных
# =============================================================================
if [[ "$OFFLINE" == false ]]; then
    echo -e "\n📥 Загрузка данных..."

    # Геном B. subtilis 168 (~4 МБ)
    if [[ ! -f "$ASSEMBLY" ]]; then
        echo "   Геном B. subtilis..."
        wget -qc --show-progress \
            -O "${ASSEMBLY}.gz" \
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz"
        gunzip -f "${ASSEMBLY}.gz"
    else
        echo "   Геном: используем кэш"
    fi

    # Swiss-Prot (первые 100k записей — ~40 МБ вместо 600 МБ)
    if [[ ! -f "$PROTEIN_DB_FASTA" ]]; then
        echo "   Swiss-Prot (подмножество)..."
        wget -qc --show-progress \
            -O "${DB_DIR}/uniprot_sprot.fasta.gz" \
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
        zcat "${DB_DIR}/uniprot_sprot.fasta.gz" \
            | awk '/^>/{n++} n>100000{exit} {print}' \
            > "$PROTEIN_DB_FASTA"
        rm "${DB_DIR}/uniprot_sprot.fasta.gz"
    else
        echo "   Swiss-Prot: используем кэш"
    fi

    # Pfam — только 10 «учебных» доменов (ABC-транспортёры, шапероны и др.)
    if [[ ! -f "${HMM_DB}.h3i" ]]; then
        echo "   Pfam-A (подмножество доменов)..."
        wget -qc --show-progress \
            -O "${DB_DIR}/Pfam-A.hmm.gz" \
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
        gunzip -kf "${DB_DIR}/Pfam-A.hmm.gz"

        # Извлекаем конкретные домены
        DOMAINS="PF00005 PF00009 PF00012 PF00013 PF00023 PF00027 PF00028 PF00043 PF00044 PF00072"
        python3 -c '
import sys
domains = set(sys.argv[1].split())
capture = False
buffer = []
with open(sys.argv[2], "r", encoding="utf-8") as f:
    for line in f:
        if line.startswith("HMMER"):
            buffer = [line]
            continue
        buffer.append(line)
        if line.startswith("ACC "):
            acc = line.split()[1].split(".")[0]
            if acc in domains:
                capture = True
        elif line.startswith("//"):
            if capture:
                sys.stdout.write("".join(buffer))
                capture = False
            buffer = []
' "$DOMAINS" "${DB_DIR}/Pfam-A.hmm" > "$HMM_DB"
        hmmpress -f "$HMM_DB"
        rm "${DB_DIR}/Pfam-A.hmm" "${DB_DIR}/Pfam-A.hmm.gz" 2>/dev/null || true
    else
        echo "   Pfam: используем кэш"
    fi
else
    echo -e "\n⚠️  Офлайн-режим: пропуск загрузки"
    for f in "$ASSEMBLY" "$PROTEIN_DB_FASTA" "${HMM_DB}.h3i"; do
        [[ -f "$f" ]] || { echo "❌ Кэш не найден: $f"; exit 1; }
    done
fi

# =============================================================================
# 3. DIAMOND — база данных
# =============================================================================
if [[ ! -f "${DIAMOND_DB}.dmnd" ]]; then
    echo -e "\n🔨 Создание DIAMOND БД..."
    diamond makedb --in "$PROTEIN_DB_FASTA" --db "$DIAMOND_DB" --threads "$CPUS" --quiet
else
    echo -e "\n🔨 DIAMOND БД: используем кэш"
fi

# =============================================================================
# 4. Prokka — структурная аннотация
# =============================================================================
echo -e "\n🧬 Prokka..."
prokka \
    --outdir "${OUT_DIR}/prokka" \
    --prefix "$PROKKA_PREFIX" \
    --kingdom Bacteria \
    --cpus "$CPUS" \
    --force \
    --quiet \
    "$ASSEMBLY"

FAA="${OUT_DIR}/prokka/${PROKKA_PREFIX}.faa"

# =============================================================================
# 5. DIAMOND blastp — функциональная гомология
# =============================================================================
echo -e "\n🔎 DIAMOND blastp..."
diamond blastp \
    --db "$DIAMOND_DB" \
    --query "$FAA" \
    --out "${OUT_DIR}/diamond/blastp.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --evalue 1e-5 \
    --max-target-seqs 3 \
    --threads "$CPUS" \
    --block-size 1.0 \
    --very-sensitive \
    --quiet

# =============================================================================
# 6. hmmscan — доменная аннотация
# =============================================================================
echo -e "\n🔬 hmmscan..."
QUERY_FAA="$FAA"
if [[ "$MAX_PROTEINS" -gt 0 ]]; then
    seqkit head -n "$MAX_PROTEINS" "$FAA" > "${OUT_DIR}/hmmer/query.faa"
    QUERY_FAA="${OUT_DIR}/hmmer/query.faa"
fi
hmmscan \
    --cpu "$CPUS" \
    --domtblout "${OUT_DIR}/hmmer/domtblout.txt" \
    --noali \
    -E 1e-5 --domE 1e-3 \
    -o "${OUT_DIR}/hmmer/hmmscan.log" \
    "$HMM_DB" \
    "$QUERY_FAA"

# =============================================================================
# 7. Итоговая статистика
# =============================================================================
echo -e "\n📊 Статистика:"
if [[ -f "${OUT_DIR}/prokka/${PROKKA_PREFIX}.txt" ]]; then
    echo "   [Prokka]"
    grep -E "CDS|tRNA|rRNA" "${OUT_DIR}/prokka/${PROKKA_PREFIX}.txt" \
        | awk '{print "     "$0}'
fi

TOTAL=$(grep -c "^>" "$FAA" || echo 0)
if [[ -f "${OUT_DIR}/diamond/blastp.tsv" ]]; then
    HITS=$(cut -f1 "${OUT_DIR}/diamond/blastp.tsv" | sort -u | wc -l | tr -d ' ')
    echo "   [DIAMOND]  $HITS / $TOTAL белков нашли гомологию"
fi

if [[ -f "${OUT_DIR}/hmmer/domtblout.txt" ]]; then
    NDOM=$(grep -cv "^#" "${OUT_DIR}/hmmer/domtblout.txt" || echo 0)
    echo "   [HMMER]    Доменных совпадений: $NDOM"
fi

echo -e "\n✅ Готово! Результаты: ${OUT_DIR}/"
