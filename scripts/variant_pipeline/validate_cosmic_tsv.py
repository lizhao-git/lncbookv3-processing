#!/usr/bin/env python3
import argparse
import csv
import shutil
import sys

REQUIRED_GROUPS = {
    "chrom": {"chrom", "chr", "chromosome"},
    "pos": {"pos", "position", "genome_position", "start"},
    "ref": {"ref", "reference", "ref_allele"},
    "alt": {"alt", "alternate", "alt_allele", "mutation_allele"},
    "cosmic_id": {"cosmic_id", "mutation_id", "mutation id", "mutationid"},
    "fathmm_mkl_score": {"fathmm_mkl_score", "fathmm-mkl_score", "fathmm score", "fathmm_mkl"},
    "disease_name": {"disease_name", "primary_histology", "histology", "tumour_type", "disease"},
}


def normalize(name: str) -> str:
    return name.strip().lower().replace("-", "_")


def main():
    parser = argparse.ArgumentParser(description="Validate COSMIC tabular variant file")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    errors = []
    record_count = 0

    with open(args.input_tsv, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            errors.append((0, "Input file is empty"))
            header = []

        normalized = {normalize(col): idx for idx, col in enumerate(header)}
        resolved = {}
        for canonical, aliases in REQUIRED_GROUPS.items():
            idx = None
            for alias in aliases:
                if alias in normalized:
                    idx = normalized[alias]
                    break
            if idx is None:
                errors.append((0, f"Missing required column for {canonical}"))
            else:
                resolved[canonical] = idx

        for ln, row in enumerate(reader, start=2):
            if not row or all(not cell.strip() for cell in row):
                continue
            record_count += 1
            if len(row) < len(header):
                errors.append((ln, "Row has fewer columns than header"))
                continue
            if resolved:
                try:
                    pos_val = row[resolved["pos"]].strip()
                    int(pos_val)
                except Exception:
                    errors.append((ln, "Position is not an integer"))
                score_val = row[resolved["fathmm_mkl_score"]].strip()
                if score_val not in {"", "NA", "."}:
                    try:
                        float(score_val)
                    except ValueError:
                        errors.append((ln, "FATHMM-MKL score is not numeric"))

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"records\t{record_count}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        if errors:
            rep.write("errors\t" + " | ".join([f"line {ln}: {msg}" for ln, msg in errors[:100]]) + "\n")

    if errors:
        sys.stderr.write("COSMIC TSV validation failed. See report for details.\n")
        sys.exit(1)

    shutil.copyfile(args.input_tsv, args.output_tsv)


if __name__ == "__main__":
    main()
