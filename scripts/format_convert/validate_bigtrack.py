#!/usr/bin/env python3
import argparse

from genomic_intervals.bigtrack import validate_bigtrack


def main():
    parser = argparse.ArgumentParser(description="Validate a bigWig/bigBed file with kent bigWigInfo/bigBedInfo")
    parser.add_argument("--input-track", required=True)
    parser.add_argument("--output-track", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_bigtrack(args.input_track, args.output_track, args.report)


if __name__ == "__main__":
    main()
