#!/usr/bin/env python3
import argparse

from genomic_intervals.annotation import validate_annotation


def main():
    parser = argparse.ArgumentParser(description="Validate GTF annotation and copy it through")
    parser.add_argument("--input-gtf", required=True)
    parser.add_argument("--output-gtf", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_annotation(args.input_gtf, args.output_gtf, args.report, fmt="gtf")


if __name__ == "__main__":
    main()
