# Genome Annotation

Functional annotation of bacterial genomes in a single script.

## Stack

| Tool | Task |
|---|---|
| **Prokka** | Structural annotation (ORF, tRNA, rRNA) |
| **DIAMOND blastp** | Functional homology mapping (Swiss-Prot) |
| **hmmscan** | Domain annotation (Pfam) |

## Usage

```bash
conda env create -f environment.yml
conda activate annotation

# Full run (downloads a test genome B. subtilis ~4 MB)
bash run_annotation.sh

# Offline mode (data already present in data/)
bash run_annotation.sh --offline

# Limit the number of proteins for hmmscan (for testing)
bash run_annotation.sh --max-proteins 200 --cpus 4
```

## Output Structure

```
outputs/run_YYYYMMDD_HHMMSS/
├── prokka/             # .gff, .faa, .ffn, .txt
├── diamond/            # blastp.tsv
└── hmmer/              # domtblout.txt
```

## Running Your Own Genome

Place your `assembly.fasta` into `data/input/` and run the script with the `--offline` flag.
