#!/usr/bin/env python3
import argparse

from genomic_intervals.vcf import validate_vcf


def main():
    parser = argparse.ArgumentParser(description="Validate a VCF file and copy it through")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_vcf(args.input_vcf, args.output_vcf, args.report)


if __name__ == "__main__":
    main()
