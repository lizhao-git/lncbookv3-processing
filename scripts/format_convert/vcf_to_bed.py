#!/usr/bin/env python3
"""Convert VCF records to BED intervals.

Each record becomes one BED row spanning the reference allele
(0-based half-open): start = POS - 1, end = POS - 1 + len(REF). The name
column holds the record ID (falling back to ``chrom_pos`` when IDs are
missing) and the QUAL value is copied into the score column. The ALT
text is appended as a seventh column.
"""
import argparse

from genomic_intervals.io import open_text
from genomic_intervals.vcf import iter_vcf


def convert_vcf_to_bed(input_path: str, output_path: str, name_field: str) -> int:
    records = 0
    with open(output_path, "w", encoding="utf-8") as out:
        for _line_no, kind, item in iter_vcf(input_path):
            if kind != "record":
                continue
            if name_field == "chrom_pos" or item.id in ("", "."):
                name = f"{item.chrom}_{item.pos}"
            else:
                name = item.id
            start = item.pos - 1
            end = start + len(item.ref)
            out.write("\t".join([
                item.chrom, str(start), str(end), name, item.qual, ".", item.alt_text(),
            ]) + "\n")
            records += 1
    return records


def main():
    parser = argparse.ArgumentParser(description="Convert VCF records to BED intervals")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--name-field", choices=["id", "chrom_pos"], default="id")
    args = parser.parse_args()
    convert_vcf_to_bed(args.input_vcf, args.output_bed, args.name_field)


if __name__ == "__main__":
    main()
