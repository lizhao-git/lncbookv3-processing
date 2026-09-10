#!/usr/bin/env python3
import argparse

from genomic_intervals.bed import validate_bed


def main():
    parser = argparse.ArgumentParser(description="Validate BED coordinate file and copy it through")
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--min-columns", type=int, default=3)
    args = parser.parse_args()
    validate_bed(args.input_bed, args.output_bed, args.report, min_columns=args.min_columns)


if __name__ == "__main__":
    main()
