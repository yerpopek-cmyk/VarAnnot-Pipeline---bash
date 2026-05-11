# Genome Annotation

Функциональная аннотация бактериального генома за один скрипт.

## Стек

| Инструмент | Задача |
|---|---|
| **Prokka** | Структурная аннотация (ORF, тРНК, рРНК) |
| **DIAMOND blastp** | Функциональная гомология (Swiss-Prot) |
| **hmmscan** | Доменная аннотация (Pfam) |

## Запуск

```bash
conda env create -f environment.yml
conda activate annotation

# Полный запуск (скачает тестовый геном B. subtilis ~4 МБ)
bash run_annotation.sh

# Офлайн (данные уже в data/)
bash run_annotation.sh --offline

# Ограничить число белков для hmmscan (тест)
bash run_annotation.sh --max-proteins 200 --cpus 4
```

## Структура выходных данных

```
outputs/run_YYYYMMDD_HHMMSS/
├── prokka/             # .gff, .faa, .ffn, .txt
├── diamond/            # blastp.tsv
└── hmmer/              # domtblout.txt
```

## Размещение своего генома

Положи `assembly.fasta` в `data/input/` и запусти с `--offline`.
