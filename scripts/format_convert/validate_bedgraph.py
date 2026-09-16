#!/usr/bin/env python3
import argparse

from genomic_intervals.bedgraph import validate_bedgraph


def main():
    parser = argparse.ArgumentParser(description="Validate a bedGraph coverage file and copy it through")
    parser.add_argument("--input-bedgraph", required=True)
    parser.add_argument("--output-bedgraph", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate_bedgraph(args.input_bedgraph, args.output_bedgraph, args.report)


if __name__ == "__main__":
    main()
