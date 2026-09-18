#!/usr/bin/env python3
import argparse

from genomic_intervals.plink import validate_plink


def main():
    parser = argparse.ArgumentParser(
        description="Validate a PLINK .map/.ped pair and copy both files through")
    parser.add_argument("--input-map", required=True)
    parser.add_argument("--input-ped", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--output-ped", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_plink(args.input_map, args.input_ped,
                   args.output_map, args.output_ped, args.report)


if __name__ == "__main__":
    main()
