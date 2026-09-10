#!/usr/bin/env python3
import argparse

from genomic_intervals.annotation import validate_annotation


def main():
    parser = argparse.ArgumentParser(description="Validate GTF/GFF3 annotation and copy it to a normalized text file")
    parser.add_argument("--input-annotation", required=True)
    parser.add_argument("--output-annotation", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--format", choices=["auto", "gtf", "gff3"], default="auto")
    args = parser.parse_args()
    validate_annotation(args.input_annotation, args.output_annotation, args.report, fmt=args.format)


if __name__ == "__main__":
    main()

