#!/usr/bin/env python3
"""
05_report.py — Приоритизация вариантов и Markdown-отчёт.

Парсит поле CSQ (Ensembl VEP) из аннотированного VCF,
вычисляет составной балл, выводит TSV и отчёт.

Использование:
    python3 05_report.py <annotated.vcf.gz> <output_dir> [--top N]
"""
import argparse
import gzip
import sys
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# Весовые таблицы
# ---------------------------------------------------------------------------
IMPACT_SCORE: dict[str, int] = {
    "HIGH": 4, "MODERATE": 2, "LOW": 1, "MODIFIER": 0
}

# ClinVar-значения приходят в нижнем регистре через VEP
CLNSIG_SCORE: dict[str, int] = {
    "pathogenic": 5,
    "likely_pathogenic": 3,
    "uncertain_significance": 0,
    "vus": 0,
    "likely_benign": -1,
    "benign": -2,
}

# ---------------------------------------------------------------------------
# Парсинг
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("vcf",    help="Аннотированный VCF (VEP, bgzip)")
    p.add_argument("outdir", help="Выходная директория")
    p.add_argument("--top",  type=int, default=20,
                   help="Число топ-вариантов в отчёте (по умолчанию: 20)")
    return p.parse_args()


def _csq_fields_from_header(header_lines: list[str]) -> list[str]:
    """Извлечь поля CSQ из строки ##INFO=<ID=CSQ,...> в заголовке VCF."""
    for line in header_lines:
        if "ID=CSQ" in line and "Format:" in line:
            fragment = line.split("Format:")[1].rstrip('">')
            return [f.strip() for f in fragment.split("|")]
    return []


def _info_value(info: str, key: str) -> str | None:
    """Вернуть значение поля KEY= из строки INFO."""
    prefix = f"{key}="
    for field in info.split(";"):
        if field.startswith(prefix):
            return field[len(prefix):]
    return None


def _compute_score(impact: str, clnsig_raw: str, af_raw: str) -> tuple[int, float | None]:
    """Вычислить приоритетный балл варианта."""
    score = IMPACT_SCORE.get(impact, 0)

    # ClinVar (несколько значений через &)
    if clnsig_raw:
        for sig in clnsig_raw.lower().replace(" ", "_").split("&"):
            score += CLNSIG_SCORE.get(sig.strip(), 0)

    # Частота gnomAD
    af: float | None = None
    if af_raw and af_raw not in ("", "."):
        try:
            af = float(af_raw.split("&")[0])
        except ValueError:
            pass
    if af is not None:
        if af < 0.001:
            score += 3
        elif af < 0.01:
            score += 1
        elif af > 0.05:
            score -= 1

    return score, af


def parse_vcf(vcf_path: str) -> list[dict]:
    """
    Распарсить аннотированный VCF.
    Возвращает список вариантов (только PASS или неразмеченные).
    """
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    header_lines: list[str] = []
    csq_fields: list[str] = []
    variants: list[dict] = []

    with open_fn(vcf_path, "rt") as fh:
        for line in fh:
            line = line.rstrip()

            if line.startswith("##"):
                header_lines.append(line)
                # Читаем поля CSQ на лету, чтобы не перечитывать заголовок
                if not csq_fields and "ID=CSQ" in line:
                    csq_fields = _csq_fields_from_header([line])
                continue

            if line.startswith("#"):
                continue  # строка с именами колонок

            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, _, ref, alt, qual, filt, info = cols[:8]

            # Пропускаем варианты с мягкими фильтрами (LowQual / LowAB)
            if filt not in (".", "PASS", ""):
                continue

            csq_raw = _info_value(info, "CSQ")
            if not csq_raw:
                continue

            # Берём первую запись CSQ — VEP сортирует по убыванию значимости
            first = csq_raw.split(",")[0].split("|")
            csq = dict(zip(csq_fields, first)) if csq_fields else {}

            impact = csq.get("IMPACT", "MODIFIER")
            gene   = csq.get("SYMBOL", ".")
            conseq = csq.get("Consequence", ".").replace("_variant", "")
            hgvsp  = csq.get("HGVSp", ".")
            af_raw = (
                csq.get("gnomADg_AF")
                or csq.get("gnomADe_AF")
                or csq.get("AF", "")
            )
            clnsig = csq.get("CLIN_SIG", "")

            score, af = _compute_score(impact, clnsig, af_raw)

            variants.append({
                "CHROM":       chrom,
                "POS":         pos,
                "REF":         ref,
                "ALT":         alt,
                "GENE":        gene,
                "CONSEQUENCE": conseq,
                "IMPACT":      impact,
                "HGVSp":       hgvsp,
                "AF":          f"{af:.5f}" if af is not None else ".",
                "CLIN_SIG":    clnsig or ".",
                "SCORE":       score,
            })

    variants.sort(key=lambda v: v["SCORE"], reverse=True)
    return variants


# ---------------------------------------------------------------------------
# Запись результатов
# ---------------------------------------------------------------------------
def write_tsv(variants: list[dict], path: Path) -> None:
    if not variants:
        path.write_text("")
        return
    keys = list(variants[0].keys())
    lines = ["\t".join(keys)]
    for v in variants:
        lines.append("\t".join(str(v[k]) for k in keys))
    path.write_text("\n".join(lines) + "\n")


def write_report(variants: list[dict], path: Path, top_n: int) -> None:
    now  = datetime.now().strftime("%Y-%m-%d %H:%M")
    top  = variants[:top_n]

    high = sum(1 for v in variants if v["IMPACT"] == "HIGH")
    mod  = sum(1 for v in variants if v["IMPACT"] == "MODERATE")
    path_var = sum(1 for v in variants if "pathogenic" in v["CLIN_SIG"].lower())

    header = "| # | GENE | CHROM:POS | REF→ALT | CONSEQUENCE | IMPACT | AF | CLIN_SIG | SCORE |"
    sep    = "|---|------|-----------|---------|-------------|--------|----|----------|-------|"
    rows = []
    for i, v in enumerate(top, 1):
        rows.append(
            f"| {i} | **{v['GENE']}** | {v['CHROM']}:{v['POS']} "
            f"| {v['REF']}→{v['ALT']} | {v['CONSEQUENCE']} "
            f"| {v['IMPACT']} | {v['AF']} | {v['CLIN_SIG']} | {v['SCORE']} |"
        )

    path.write_text(f"""\
# Отчёт по вариантам
Создан: {now}

## Сводка

| Показатель | Значение |
|---|---|
| Всего PASS-вариантов | {len(variants)} |
| HIGH impact | {high} |
| MODERATE impact | {mod} |
| Патогенных (ClinVar) | {path_var} |

## Топ-{top_n} по приоритету

{header}
{sep}
{chr(10).join(rows)}

## Формула балла

```
SCORE = IMPACT_score + gnomAD_AF_score + ClinVar_score

IMPACT:   HIGH=4, MODERATE=2, LOW=1, MODIFIER=0
gnomAD AF: <0.001 → +3  |  0.001–0.01 → +1  |  >0.05 → −1
ClinVar:  Pathogenic=+5, Likely_pathogenic=+3, Uncertain=0,
          Likely_benign=−1, Benign=−2
```
""")


# ---------------------------------------------------------------------------
# Точка входа
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sample = Path(args.vcf).name.split(".")[0]

    print(f"[05] Парсинг: {args.vcf}")
    variants = parse_vcf(args.vcf)
    print(f"[05] PASS-вариантов: {len(variants)}")

    if not variants:
        print("[05] ⚠️  Нет вариантов для отчёта.", file=sys.stderr)
        return

    tsv_path = outdir / f"{sample}.prioritized.tsv"
    md_path  = outdir / f"{sample}.report.md"

    write_tsv(variants, tsv_path)
    write_report(variants, md_path, args.top)

    print(f"[05] TSV    : {tsv_path}")
    print(f"[05] Отчёт : {md_path}")


if __name__ == "__main__":
    main()
