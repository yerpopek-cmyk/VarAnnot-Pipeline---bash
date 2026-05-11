# Variant Pipeline

Полный пайплайн от FASTQ до приоритизированного отчёта по вариантам.

## Стек

| Шаг | Инструмент | Задача |
|-----|-----------|--------|
| 01  | BWA-MEM2 + samtools | Выравнивание, сортировка, MarkDup |
| 02  | FreeBayes | SNV/indel calling |
| 03  | bcftools  | Нормализация + soft-фильтрация |
| 04  | Ensembl VEP | Аннотация: последствия + gnomAD AF + ClinVar |
| 05  | Python 3  | Приоритизация + Markdown-отчёт |

> **Почему VEP вместо snpEff + отдельных bcftools annotate?**  
> VEP за один запуск даёт функциональные последствия, gnomAD AF и ClinVar статус
> — это три шага в одном, без отдельных скачиваний gnomAD VCF (~50 ГБ) и ClinVar.

## Быстрый старт

```bash
# 1. Окружение
conda env create -f env/environment.yml
conda activate variants

# 2. Кэш VEP (один раз, ~15 ГБ)
vep_install -a cf -s homo_sapiens -y GRCh38 \
            --CACHEDIR data/db/vep_cache --NO_HTSLIB

# 3. Настроить пути
cp workflow/config.sh.example workflow/config.sh
nano workflow/config.sh   # указать REF_FASTA, READS_R1/R2, SAMPLE_ID

# 4. Запуск
bash workflow/run_all.sh
```

## Структура выходных данных

```
outputs/run_YYYYMMDD_HHMMSS/
├── 1_bams/               # sorted + markdup BAM
├── 2_vcf_raw/            # сырой VCF (FreeBayes)
├── 3_vcf_filtered/       # нормализованный, soft-filtered, аннотированный VCF
└── 4_reports/
    ├── sample.prioritized.tsv   # все PASS-варианты с баллами
    └── sample.report.md         # топ-N с формулой балла
```

## Перезапуск с определённого шага

```bash
# Пересчитать только отчёт (шаг 05)
RUN_DIR=outputs/run_20240501_120000 bash workflow/run_all.sh --from 05
```

## Если у тебя уже есть BAM

Укажи в `config.sh`:
```bash
BAM_INPUT="data/input/sample.bam"
```
и запусти с `--from 02`.
