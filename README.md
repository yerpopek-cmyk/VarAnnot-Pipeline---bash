# NGS Pipeline: Comprehensive Genomics Workflow

![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-blue)
![Bash](https://img.shields.io/badge/Language-Bash%20%7C%20Python-green)
![License](https://img.shields.io/badge/License-MIT-orange)

Welcome to the **NGS Pipeline** repository! This project provides two robust, independent, and automated bioinformatics pipelines designed for genomic data processing. Whether you are performing functional annotation of bacterial genomes or calling clinically relevant variants in human data, this repository offers streamlined, reproducible workflows.

## 🧬 Pipeline Architectures

### 1. `annotation/` — Bacterial Genome Annotation
A comprehensive functional annotation workflow optimized for speed and accuracy.
- **Prokka**: Performs rapid structural annotation (ORFs, tRNAs, rRNAs).
- **DIAMOND blastp**: Fast functional homology mapping against Swiss-Prot.
- **HMMER (hmmscan)**: Deep domain annotation using the Pfam database.

### 2. `variants/` — Human Variant Calling & Interpretation
A clinically-oriented pipeline transforming raw FASTQ reads into annotated, prioritized variants.
- **BWA-MEM2 & Samtools**: Read alignment, sorting, and PCR duplicate marking.
- **FreeBayes**: Highly accurate haplotype-based variant calling.
- **BCFtools**: VCF normalization, splitting multiallelic sites, and soft-filtering.
- **Ensembl VEP**: Exhaustive functional annotation (Consequences, gnomAD AF, ClinVar significance).
- **Custom Python Reporter**: Prioritizes variants based on clinical impact and allele balance.

---

## 🚀 Quick Start Guide

### Prerequisites
- **OS**: Linux or macOS (Windows Subsystem for Linux / WSL2 supported)
- **Environment Manager**: `Conda` or `Micromamba`
- **Storage**: ~20 GB free space (primarily for the VEP cache and reference databases)

### Option A: Bacterial Genome Annotation

```bash
# 1. Navigate to the annotation directory
cd annotation

# 2. Create and activate the conda environment
conda env create -f environment.yml
conda activate annotation

# 3. Run the pipeline (downloads test data automatically)
bash run_annotation.sh
```

### Option B: Human Variant Calling Pipeline

```bash
# 1. Navigate to the variants directory
cd variants

# 2. Create and activate the environment
conda env create -f env/environment.yml
conda activate variants

# 3. Setup configuration
cp workflow/config.sh.example workflow/config.sh
# Note: Edit workflow/config.sh to set your FASTQ paths and sample IDs!

# 4. Execute the full orchestrator
bash workflow/run_all.sh
```

---

## 🛠️ Repository Structure

```text
ngs_pipeline/
├── annotation/             # Bacterial annotation pipeline (Prokka, DIAMOND, HMMER)
│   ├── data/               # Databases and input/output directories
│   ├── environment.yml     # Conda dependencies
│   └── run_annotation.sh   # Main execution script
├── variants/               # Human variant calling pipeline (BWA, FreeBayes, VEP)
│   ├── env/                # Conda dependencies
│   ├── scripts/            # Modular steps (01_align, 02_call, etc.)
│   └── workflow/           # Orchestrator and configuration
├── .gitignore              # Pre-configured to ignore massive genomics data
├── Makefile                # Quick alias commands for standard runs
└── README.md               # This documentation file
```

## 📊 Output Interpretations

Each pipeline generates neatly organized output directories timestamped for reproducibility:
- **Annotation Outputs**: Found in `annotation/outputs/run_YYYYMMDD_HHMMSS/`, containing standard GFF3, FASTA, and tabular homology reports.
- **Variant Outputs**: Found in `variants/outputs/run_YYYYMMDD_HHMMSS/`, containing intermediate BAMs, raw/filtered VCFs, and the final clinically prioritized CSV report.

## 🤝 Contributing
Contributions, issues, and feature requests are welcome! Feel free to check the issues page.

---
*Built with ❤️ for the bioinformatics community.*
