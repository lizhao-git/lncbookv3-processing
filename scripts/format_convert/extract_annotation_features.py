#!/usr/bin/env python3
import argparse

from genomic_intervals.annotation import extract_features


def main():
    parser = argparse.ArgumentParser(description="Extract annotation features from GTF/GFF3 into a BED-like table")
    parser.add_argument("--input-annotation", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--format", choices=["auto", "gtf", "gff3"], default="auto")
    parser.add_argument("--feature", nargs="+", choices=["gene", "transcript", "exon", "intron"])
    args = parser.parse_args()
    extract_features(args.input_annotation, args.output_bed, fmt=args.format, features=args.feature)


if __name__ == "__main__":
    main()
