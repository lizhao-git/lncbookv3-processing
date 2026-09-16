#!/usr/bin/env python3
import argparse

from genomic_intervals.liftover import run_liftover, write_report


def main():
    parser = argparse.ArgumentParser(description="Lift BED-family intervals across assemblies with kent liftOver")
    parser.add_argument("--input-track", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--unmapped-bed", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--min-match", type=float, default=0.95)
    args = parser.parse_args()
    metrics = run_liftover(args.input_track, args.chain, args.output_bed, args.unmapped_bed, min_match=args.min_match)
    write_report(args.report, args.input_track, args.chain, args.min_match, metrics)


if __name__ == "__main__":
    main()
