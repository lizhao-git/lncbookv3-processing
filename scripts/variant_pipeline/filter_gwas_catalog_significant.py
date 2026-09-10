#!/usr/bin/env python3
import argparse
import csv
from collections import Counter

PVALUE_THRESHOLD = 5e-8
COLUMN_ALIASES = {
    "chrom": ["chrom", "chr", "chromosome", "chr_id"],
    "pos": ["pos", "position", "base_pair_location", "bp", "start", "chr_pos"],
    "ref": ["ref", "reference", "reference_allele", "other_allele"],
    "alt": ["alt", "alternate", "alt_allele", "risk_allele", "effect_allele", "strongest_snp_risk_allele"],
    "variant_id": ["variant_id", "rsid", "snp", "snps", "markername", "snp_id_current"],
    "p_value": ["p_value", "p-value", "pvalue", "p value"],
    "trait_name": ["trait_name", "disease_trait", "trait", "mapped_trait", "disease/trait"],
    "mapped_trait": ["mapped_trait"],
    "mapped_trait_uri": ["mapped_trait_uri"],
    "risk_allele_frequency": ["risk_allele_frequency"],
    "or_or_beta": ["or_or_beta", "or or beta"],
    "ci_95_text": ["95%_ci_(text)", "95% ci (text)", "ci", "confidence_interval", "95_ci_text"],
}


def normalize(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_")


def resolve_columns(fieldnames):
    normalized = {normalize(col): col for col in fieldnames}
    resolved = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            alias_norm = normalize(alias)
            if alias_norm in normalized:
                resolved[canonical] = normalized[alias_norm]
                break
    return resolved


def normalize_chrom(chrom: str) -> str:
    value = chrom.strip()
    if value.startswith("chr"):
        return value
    if value == "MT":
        return "chrM"
    return f"chr{value}"


def extract_alt_allele(raw_value: str) -> str:
    value = raw_value.strip()
    if not value:
        return "."
    if "-" in value:
        candidate = value.rsplit("-", 1)[-1].strip()
        if candidate:
            return candidate
    return value


def first_token(value: str) -> str:
    if not value:
        return "NA"
    return value.split(",", 1)[0].strip() or "NA"


def extract_efo_id(mapped_trait_uri: str) -> str:
    token = first_token(mapped_trait_uri)
    if token == "NA":
        return "NA"
    tail = token.rstrip("/").rsplit("/", 1)[-1]
    return tail or "NA"


def main():
    parser = argparse.ArgumentParser(description="Filter GWAS Catalog associations with p-value < 5e-8")
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
        required = ["chrom", "pos", "variant_id", "p_value", "trait_name"]
        missing = [name for name in required if name not in resolved]
        if missing:
            raise SystemExit(f"Missing required GWAS Catalog columns: {', '.join(missing)}")

        with open(args.output_bed, "w", encoding="utf-8") as bed, open(args.output_tsv, "w", encoding="utf-8", newline="") as out_tsv:
            writer = csv.DictWriter(out_tsv, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()

            for row in reader:
                p_value_text = str(row[resolved["p_value"]]).strip()
                try:
                    p_value = float(p_value_text)
                except ValueError:
                    continue
                if p_value >= PVALUE_THRESHOLD:
                    continue

                chrom_raw = str(row[resolved["chrom"]]).strip()
                pos_raw = str(row[resolved["pos"]]).strip()
                if not chrom_raw or not pos_raw:
                    continue
                try:
                    pos = int(pos_raw)
                except ValueError:
                    continue
                chrom = normalize_chrom(chrom_raw)
                ref = str(row[resolved["ref"]]).strip() if "ref" in resolved else "N"
                raw_alt = str(row[resolved["alt"]]).strip() if "alt" in resolved else "."
                alt = extract_alt_allele(raw_alt)
                variant_id = str(row[resolved["variant_id"]]).strip() or f"{chrom}:{pos}:{ref}:{alt}"
                trait_name = str(row[resolved["trait_name"]]).strip() or "NA"
                mapped_trait = str(row[resolved["mapped_trait"]]).strip() if "mapped_trait" in resolved else "NA"
                mapped_trait = first_token(mapped_trait)
                mapped_trait_uri = str(row[resolved["mapped_trait_uri"]]).strip() if "mapped_trait_uri" in resolved else "NA"
                efo_id = extract_efo_id(mapped_trait_uri)
                risk_allele_frequency = str(row[resolved["risk_allele_frequency"]]).strip() if "risk_allele_frequency" in resolved else "NA"
                risk_allele_frequency = risk_allele_frequency or "NA"
                or_or_beta = str(row[resolved["or_or_beta"]]).strip() if "or_or_beta" in resolved else "NA"
                or_or_beta = or_or_beta or "NA"
                ci_95_text = str(row[resolved["ci_95_text"]]).strip() if "ci_95_text" in resolved else "NA"
                ci_95_text = ci_95_text or "NA"
                start0 = pos - 1
                end = start0 + max(1, len(ref))
                label = "Genome-wide significant"

                bed.write(
                    f"{chrom}\t{start0}\t{end}\t{variant_id}\t{label}\t{variant_id}\t{trait_name}\t{p_value}\tNA\t"
                    f"{mapped_trait}\t{efo_id}\t{risk_allele_frequency}\t{or_or_beta}\t{ci_95_text}\n"
                )
                writer.writerow(row)
                kept += 1
                counts[label] += 1

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"kept_associations\t{kept}\n")
        rep.write(f"pvalue_threshold\t{PVALUE_THRESHOLD}\n")
        rep.write(f"label_Genome-wide significant\t{counts.get('Genome-wide significant', 0)}\n")


if __name__ == "__main__":
    main()
