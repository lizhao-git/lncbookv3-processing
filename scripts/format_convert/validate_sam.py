#!/usr/bin/env python3
import argparse

from genomic_intervals.sam import validate_sam


def main():
    parser = argparse.ArgumentParser(description="Validate a SAM file and copy it through")
    parser.add_argument("--input-sam", required=True)
    parser.add_argument("--output-sam", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_sam(args.input_sam, args.output_sam, args.report)


if __name__ == "__main__":
    main()
