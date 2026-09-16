#!/usr/bin/env python3
import argparse

from genomic_intervals.wig import validate_wig


def main():
    parser = argparse.ArgumentParser(description="Validate a WIG coverage file and copy it through")
    parser.add_argument("--input-wig", required=True)
    parser.add_argument("--output-wig", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_wig(args.input_wig, args.output_wig, args.report)


if __name__ == "__main__":
    main()
