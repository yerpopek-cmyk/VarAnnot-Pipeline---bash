# NGS Pipeline

Два независимых пайплайна для работы с геномными данными.

```
ngs_pipeline/
├── annotation/   — аннотация бактериальных геномов (Prokka + DIAMOND + HMMER)
├── variants/     — вызов и интерпретация вариантов человека (BWA → FreeBayes → VEP)
└── Makefile      — быстрые команды для обоих пайплайнов
```

## Быстрый старт

```bash
# Аннотация генома
cd annotation && conda env create -f environment.yml && conda activate annotation
bash run_annotation.sh

# Вариантный пайплайн
cd variants && conda env create -f env/environment.yml && conda activate variants
cp workflow/config.sh.example workflow/config.sh  # отредактируй пути
bash workflow/run_all.sh
```

## Требования
- Linux / macOS (WSL2)
- [Micromamba](https://mamba.readthedocs.io/) или Conda
- ~20 ГБ свободного места (VEP-кэш GRCh38 + БД)
