#!/usr/bin/env python3
import argparse
import os
import sys

if __package__ is None or __package__ == "":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from conservation.common import parse_number, read_tsv, write_tsv


def as_bool(value):
    return "1" if value else "0"


def classify_rows(rows, species, min_mapped_length, min_transcript_coverage, min_intron_coverage):
    out = []
    for row in rows:
        mapped = parse_number(row.get("true_mapped_length") or row.get("mapped_length"))
        transcript_len = parse_number(row.get("transcript_length"))
        transcript_coverage = parse_number(row.get("transcript_coverage"))
        if not transcript_coverage and transcript_len:
            transcript_coverage = mapped / transcript_len
        intron_coverage = parse_number(row.get("intron_coverage"), None)

        pass_length = mapped >= min_mapped_length
        pass_tx = transcript_coverage >= min_transcript_coverage
        pass_intron = True if intron_coverage is None else intron_coverage >= min_intron_coverage

        record = dict(row)
        record["target_species"] = species
        record["transcript_coverage"] = transcript_coverage
        record["pass_length_filter"] = as_bool(pass_length)
        record["pass_transcript_coverage_filter"] = as_bool(pass_tx)
        record["pass_intron_coverage_filter"] = as_bool(pass_intron)
        record["is_conserved"] = as_bool(pass_length and pass_tx and pass_intron)
        out.append(record)
    return out


def main():
    ap = argparse.ArgumentParser(description="Classify conservation evidence with LncBook-style thresholds")
    ap.add_argument("--statistics", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--min-mapped-length", type=float, default=50)
    ap.add_argument("--min-transcript-coverage", type=float, default=0.2)
    ap.add_argument("--min-intron-coverage", type=float, default=0.0)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows = classify_rows(
        read_tsv(args.statistics),
        args.species,
        args.min_mapped_length,
        args.min_transcript_coverage,
        args.min_intron_coverage,
    )
    preferred = [
        "target_species", "ids", "grouped_ids", "geneID", "gene_name",
        "transcript_length", "true_mapped_length", "mapped_length",
        "matched_length", "match_ratio", "transcript_coverage",
        "exon_covered_ratio", "intron_coverage", "pass_length_filter",
        "pass_transcript_coverage_filter", "pass_intron_coverage_filter",
        "is_conserved",
    ]
    fields = preferred + sorted({key for row in rows for key in row if key not in preferred})
    write_tsv(args.output, rows, fields)


if __name__ == "__main__":
    main()
