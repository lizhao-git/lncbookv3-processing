#!/usr/bin/env python3
import argparse

from genomic_intervals.annotation import extract_features


def main():
    parser = argparse.ArgumentParser(description="Extract gene/transcript/exon/intron features from GTF into BED")
    parser.add_argument("--input-gtf", required=True)
    parser.add_argument("--output-bed", required=True)
    args = parser.parse_args()
    extract_features(args.input_gtf, args.output_bed, fmt="gtf")


if __name__ == "__main__":
    main()
