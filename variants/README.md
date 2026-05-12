# Variant Pipeline

A complete pipeline starting from FASTQ to a prioritized variant report.

## Stack

| Step | Tool | Task |
|-----|-----------|--------|
| 01  | BWA-MEM2 + samtools | Alignment, sorting, MarkDup |
| 02  | FreeBayes | SNV/indel calling |
| 03  | bcftools  | Normalization + soft-filtering |
| 04  | Ensembl VEP | Annotation: consequences + gnomAD AF + ClinVar |
| 05  | Python 3  | Prioritization + Markdown report generation |

> **Why use VEP instead of snpEff + separate bcftools annotate?**  
> VEP provides functional consequences, gnomAD AF, and ClinVar status in a single run.
> This combines three steps into one, eliminating the need to download large databases like the gnomAD VCF (~50 GB) and ClinVar separately.

## Quick Start

```bash
# 1. Environment setup
conda env create -f env/environment.yml
conda activate variants

# 2. VEP Cache (only needed once, ~15 GB)
vep_install -a cf -s homo_sapiens -y GRCh38 \
            --CACHEDIR data/db/vep_cache --NO_HTSLIB

# 3. Configure paths
cp workflow/config.sh.example workflow/config.sh
nano workflow/config.sh   # Set REF_FASTA, READS_R1/R2, and SAMPLE_ID

# 4. Run the pipeline
bash workflow/run_all.sh
```

## Output Structure

```
outputs/run_YYYYMMDD_HHMMSS/
├── 1_bams/               # sorted + markdup BAM files
├── 2_vcf_raw/            # raw VCF (FreeBayes)
├── 3_vcf_filtered/       # normalized, soft-filtered, annotated VCF
└── 4_reports/
    ├── sample.prioritized.tsv   # All PASS variants with scores
    └── sample.report.md         # Top N variants with scoring breakdown
```

## Restarting from a specific step

```bash
# Rerun only the report generation (Step 05)
RUN_DIR=outputs/run_20240501_120000 bash workflow/run_all.sh --from 05
```

## If you already have a BAM file

Define this in `config.sh`:
```bash
BAM_INPUT="data/input/sample.bam"
```
And run the pipeline starting from step 02:
`bash workflow/run_all.sh --from 02`
