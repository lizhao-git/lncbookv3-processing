#!/usr/bin/env python3
import argparse

from genomic_intervals.chain import validate_chain


def main():
    parser = argparse.ArgumentParser(description="Validate a UCSC chain file and copy it through")
    parser.add_argument("--input-chain", required=True)
    parser.add_argument("--output-chain", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_chain(args.input_chain, args.output_chain, args.report)


if __name__ == "__main__":
    main()
