#!/usr/bin/env python3
import argparse

from genomic_intervals.genepred import validate_genepred


def main():
    parser = argparse.ArgumentParser(description="Validate a genePred/refFlat table and copy it through")
    parser.add_argument("--input-genepred", required=True)
    parser.add_argument("--output-genepred", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_genepred(args.input_genepred, args.output_genepred, args.report)


if __name__ == "__main__":
    main()
