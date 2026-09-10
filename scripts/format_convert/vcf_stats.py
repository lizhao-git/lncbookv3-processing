#!/usr/bin/env python3
"""Summarize VCF content into a metric/value report.

Record counts are site-level; allele type counts are allele-level (each
ALT allele counted once). Ts/Tv considers SNV alleles only and is
reported as NA when no transversions are present. Per-contig record
counts are emitted as ``contig.<name>`` rows.
"""
import argparse
from collections import Counter

from genomic_intervals.vcf import TRANSITIONS, allele_type, iter_vcf


def summarize_vcf(input_path: str, report_path: str):
    records = 0
    pass_records = 0
    multiallelic_records = 0
    no_alt_records = 0
    alt_alleles = 0
    type_counts = Counter()
    contig_counts = Counter()
    transitions = 0
    transversions = 0

    for _line_no, kind, item in iter_vcf(input_path):
        if kind != "record":
            continue
        records += 1
        contig_counts[item.chrom] += 1
        if item.filt == "PASS":
            pass_records += 1
        if len(item.alts) > 1:
            multiallelic_records += 1
        if not item.alts:
            no_alt_records += 1
        for alt in item.alts:
            alt_alleles += 1
            kind_name = allele_type(item.ref, alt)
            type_counts[kind_name] += 1
            if kind_name == "snv":
                if (item.ref.upper(), alt.upper()) in TRANSITIONS:
                    transitions += 1
                else:
                    transversions += 1

    ts_tv = f"{transitions / transversions:.3f}" if transversions else "NA"
    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"records\t{records}\n")
        rep.write(f"records_pass\t{pass_records}\n")
        rep.write(f"records_no_alt\t{no_alt_records}\n")
        rep.write(f"records_multiallelic\t{multiallelic_records}\n")
        rep.write(f"alt_alleles\t{alt_alleles}\n")
        rep.write(f"alleles_snv\t{type_counts['snv']}\n")
        rep.write(f"alleles_mnv\t{type_counts['mnv']}\n")
        rep.write(f"alleles_indel\t{type_counts['indel']}\n")
        rep.write(f"alleles_other\t{type_counts['other']}\n")
        rep.write(f"transitions\t{transitions}\n")
        rep.write(f"transversions\t{transversions}\n")
        rep.write(f"ts_tv_ratio\t{ts_tv}\n")
        for chrom in sorted(contig_counts):
            rep.write(f"contig.{chrom}\t{contig_counts[chrom]}\n")


def main():
    parser = argparse.ArgumentParser(description="Summarize VCF content")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    summarize_vcf(args.input_vcf, args.report)


if __name__ == "__main__":
    main()
