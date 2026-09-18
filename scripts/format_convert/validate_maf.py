#!/usr/bin/env python3
import argparse

from genomic_intervals.maf import validate_maf


def main():
    parser = argparse.ArgumentParser(
        description="Validate a TCGA MAF file and copy it through")
    parser.add_argument("--input-maf", required=True)
    parser.add_argument("--output-maf", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_maf(args.input_maf, args.output_maf, args.report)


if __name__ == "__main__":
    main()
