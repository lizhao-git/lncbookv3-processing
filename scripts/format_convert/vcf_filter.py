#!/usr/bin/env python3
"""Filter VCF records by region, quality, depth, allele type, or FILTER status.

Records are evaluated in a fixed order (region, QUAL, INFO/DP, allele
type, FILTER) and dropped at the first failing check, so the per-reason
counters in the report never overlap. Region matching uses the POS field
(site convention); ``--region CHROM`` keeps a whole contig and
``--region CHROM:START-END`` uses 1-based inclusive coordinates.
"""
import argparse
import re

from genomic_intervals.io import open_text_write
from genomic_intervals.vcf import iter_vcf, parse_info

REGION_RE = re.compile(r"^([^:]+)(?::([0-9,]+)-([0-9,]+))?$")


def parse_region(text: str):
    match = REGION_RE.match(text.strip())
    if not match:
        raise SystemExit(f"Invalid region: {text!r}; expected CHROM or CHROM:START-END")
    chrom, start, end = match.groups()
    start_i = int(start.replace(",", "")) if start else None
    end_i = int(end.replace(",", "")) if end else None
    return chrom, start_i, end_i


def matches_regions(record, regions) -> bool:
    for chrom, start, end in regions:
        if record.chrom != chrom:
            continue
        if start is None or start <= record.pos <= end:
            return True
    return False


def record_depth(record):
    for key, value in parse_info(record.info):
        if key == "DP" and value is not None:
            try:
                return int(value)
            except ValueError:
                return None
    return None


def filter_vcf(input_path: str, output_vcf: str, report_path: str, regions,
               min_qual=None, min_dp=None, keep_types=None, require_pass=False):
    stats = {
        "input_records": 0,
        "kept_records": 0,
        "removed_region": 0,
        "removed_qual": 0,
        "removed_depth": 0,
        "removed_type": 0,
        "removed_filter": 0,
    }
    with open_text_write(output_vcf) as out:
        for _line_no, kind, item in iter_vcf(input_path):
            if kind in ("meta", "header"):
                out.write(item + "\n")
                continue
            stats["input_records"] += 1
            if regions and not matches_regions(item, regions):
                stats["removed_region"] += 1
                continue
            if min_qual is not None:
                try:
                    qual = float(item.qual)
                except ValueError:
                    qual = None
                if qual is None or qual < min_qual:
                    stats["removed_qual"] += 1
                    continue
            if min_dp is not None:
                depth = record_depth(item)
                if depth is None or depth < min_dp:
                    stats["removed_depth"] += 1
                    continue
            if keep_types and not (item.allele_types() & keep_types):
                stats["removed_type"] += 1
                continue
            if require_pass and item.filt != "PASS":
                stats["removed_filter"] += 1
                continue
            out.write(item.to_line() + "\n")
            stats["kept_records"] += 1

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        for key, value in stats.items():
            rep.write(f"{key}\t{value}\n")
    return stats


def main():
    parser = argparse.ArgumentParser(description="Filter VCF records")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--region", action="append", default=[],
                        help="CHROM or CHROM:START-END (1-based inclusive); repeatable")
    parser.add_argument("--min-qual", type=float, default=None)
    parser.add_argument("--min-dp", type=int, default=None, help="minimum INFO/DP")
    parser.add_argument("--keep-type", action="append", default=[],
                        choices=["snv", "mnv", "indel", "other"],
                        help="keep records with at least one allele of these types; repeatable")
    parser.add_argument("--require-pass", action="store_true",
                        help="keep only records whose FILTER column is PASS")
    args = parser.parse_args()

    regions = [parse_region(text) for text in args.region]
    keep_types = set(args.keep_type) or None
    filter_vcf(args.input_vcf, args.output_vcf, args.report, regions,
               min_qual=args.min_qual, min_dp=args.min_dp,
               keep_types=keep_types, require_pass=args.require_pass)


if __name__ == "__main__":
    main()
