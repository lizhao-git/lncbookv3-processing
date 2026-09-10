#!/usr/bin/env python3
import argparse
import csv
from collections import Counter

COLUMN_ALIASES = {
    "chrom": ["chrom", "chr", "chromosome"],
    "pos": ["pos", "position", "genome_position", "start"],
    "ref": ["ref", "reference", "ref_allele"],
    "alt": ["alt", "alternate", "alt_allele", "mutation_allele"],
    "cosmic_id": ["cosmic_id", "mutation_id", "mutation id", "mutationid"],
    "fathmm_mkl_score": ["fathmm_mkl_score", "fathmm-mkl_score", "fathmm score", "fathmm_mkl"],
    "disease_name": ["disease_name", "primary_histology", "histology", "tumour_type", "disease"],
}


def normalize(name: str) -> str:
    return name.strip().lower().replace("-", "_")


def resolve_columns(fieldnames):
    normalized = {normalize(col): col for col in fieldnames}
    resolved = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in normalized:
                resolved[canonical] = normalized[alias]
                break
    return resolved


def normalize_chrom(chrom: str) -> str:
    c = chrom.strip()
    if c.startswith("chr"):
        return c
    if c == "MT":
        return "chrM"
    return f"chr{c}"


def main():
    parser = argparse.ArgumentParser(description="Filter COSMIC variants with FATHMM-MKL score > 0.7")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    kept = 0
    counts = Counter()

    with open(args.input_tsv, "r", encoding="utf-8") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        resolved = resolve_columns(reader.fieldnames or [])
        required = ["chrom", "pos", "ref", "alt", "cosmic_id", "fathmm_mkl_score", "disease_name"]
        missing = [name for name in required if name not in resolved]
        if missing:
            raise SystemExit(f"Missing required COSMIC columns: {', '.join(missing)}")

        with open(args.output_bed, "w", encoding="utf-8") as bed, open(args.output_tsv, "w", encoding="utf-8", newline="") as out_tsv:
            writer = csv.DictWriter(out_tsv, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()

            for row in reader:
                score_text = str(row[resolved["fathmm_mkl_score"]]).strip()
                if score_text in {"", ".", "NA"}:
                    continue
                try:
                    score = float(score_text)
                except ValueError:
                    continue
                if score <= 0.7:
                    continue

                chrom = normalize_chrom(str(row[resolved["chrom"]]))
                pos = int(str(row[resolved["pos"]]))
                ref = str(row[resolved["ref"]]).strip() or "N"
                alt = str(row[resolved["alt"]]).strip() or "."
                cosmic_id = str(row[resolved["cosmic_id"]]).strip() or f"{chrom}:{pos}:{ref}:{alt}"
                disease_name = str(row[resolved["disease_name"]]).strip() or "NA"
                start0 = pos - 1
                end = start0 + max(1, len(ref))

                bed.write(f"{chrom}\t{start0}\t{end}\t{cosmic_id}\tPathogenic\t{cosmic_id}\t{disease_name}\t{score}\n")
                writer.writerow(row)
                kept += 1
                counts["Pathogenic"] += 1

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"kept_variants\t{kept}\n")
        rep.write(f"label_Pathogenic\t{counts.get('Pathogenic', 0)}\n")


if __name__ == "__main__":
    main()
