#!/usr/bin/env python3
import argparse

from genomic_intervals.fasta import validate_fasta


def main():
    parser = argparse.ArgumentParser(description="Validate a FASTA file and copy it through")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_fasta(args.input_fasta, args.output_fasta, args.report)


if __name__ == "__main__":
    main()
