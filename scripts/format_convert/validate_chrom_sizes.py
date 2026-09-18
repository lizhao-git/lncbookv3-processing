#!/usr/bin/env python3
import argparse

from genomic_intervals.chrom_sizes import validate_chrom_sizes


def main():
    parser = argparse.ArgumentParser(
        description="Validate a two-column chrom.sizes file and copy it through")
    parser.add_argument("--input-sizes", required=True)
    parser.add_argument("--output-sizes", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_chrom_sizes(args.input_sizes, args.output_sizes, args.report)


if __name__ == "__main__":
    main()
