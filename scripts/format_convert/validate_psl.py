#!/usr/bin/env python3
import argparse

from genomic_intervals.psl import validate_psl


def main():
    parser = argparse.ArgumentParser(description="Validate a PSL alignment file and copy it through")
    parser.add_argument("--input-psl", required=True)
    parser.add_argument("--output-psl", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--target-sizes", default=None,
                        help="chrom.sizes file; pslCheck then verifies target names and coordinates")
    parser.add_argument("--query-sizes", default=None,
                        help="query sizes file; pslCheck then verifies query names and coordinates")
    parser.add_argument("--skip-pslcheck", action="store_true",
                        help="run only the stdlib structural pass (no kent pslCheck)")
    args = parser.parse_args()
    validate_psl(args.input_psl, args.output_psl, args.report,
                 target_sizes=args.target_sizes, query_sizes=args.query_sizes,
                 use_pslcheck=not args.skip_pslcheck)


if __name__ == "__main__":
    main()
