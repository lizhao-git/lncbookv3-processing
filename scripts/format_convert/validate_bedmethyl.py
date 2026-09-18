#!/usr/bin/env python3
import argparse

from genomic_intervals.bedmethyl import validate_bedmethyl


def main():
    parser = argparse.ArgumentParser(
        description="Validate a bedMethyl file (modkit or ENCODE layout) and copy it through")
    parser.add_argument("--input-bedmethyl", required=True)
    parser.add_argument("--output-bedmethyl", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_bedmethyl(args.input_bedmethyl, args.output_bedmethyl, args.report)


if __name__ == "__main__":
    main()
