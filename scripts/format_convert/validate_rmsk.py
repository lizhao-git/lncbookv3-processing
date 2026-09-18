#!/usr/bin/env python3
import argparse

from genomic_intervals.rmsk import validate_rmsk


def main():
    parser = argparse.ArgumentParser(
        description="Validate a RepeatMasker .out/.cat file and copy it through")
    parser.add_argument("--input-rmsk", required=True)
    parser.add_argument("--output-rmsk", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_rmsk(args.input_rmsk, args.output_rmsk, args.report)


if __name__ == "__main__":
    main()
