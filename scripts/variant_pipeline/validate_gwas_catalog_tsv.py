#!/usr/bin/env python3
import argparse
import csv
import shutil
import sys

REQUIRED_GROUPS = {
    "chrom": {"chrom", "chr", "chromosome", "chr_id"},
    "pos": {"pos", "position", "base_pair_location", "bp", "start", "chr_pos"},
    "variant_id": {"variant_id", "rsid", "snp", "snps", "markername", "snp_id_current"},
    "risk_allele": {"alt", "alternate", "alt_allele", "risk_allele", "effect_allele", "strongest_snp_risk_allele"},
    "p_value": {"p_value", "p-value", "pvalue", "p value"},
    "trait_name": {"trait_name", "disease_trait", "trait", "mapped_trait", "disease/trait"},
}


def normalize(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_")


def main():
    parser = argparse.ArgumentParser(description="Validate GWAS Catalog tabular association file")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    errors = []
    warnings = []
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
                alias_norm = normalize(alias)
                if alias_norm in normalized:
                    idx = normalized[alias_norm]
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
                    int(row[resolved["pos"]].strip())
                except Exception:
                    warnings.append((ln, "Position is not an integer"))
                try:
                    float(row[resolved["p_value"]].strip())
                except Exception:
                    warnings.append((ln, "p-value is not numeric"))

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"records\t{record_count}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join([f"line {ln}: {msg}" for ln, msg in errors[:100]]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join([f"line {ln}: {msg}" for ln, msg in warnings[:100]]) + "\n")

    if errors:
        sys.stderr.write("GWAS Catalog TSV validation failed. See report for details.\n")
        sys.exit(1)

    shutil.copyfile(args.input_tsv, args.output_tsv)


if __name__ == "__main__":
    main()
