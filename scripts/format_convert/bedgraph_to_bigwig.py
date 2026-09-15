#!/usr/bin/env python3
"""Convert a bedGraph file into a bigWig track using the kent tool.

Rows are sorted into the chromosome order declared by ``--chrom-sizes``
(``track``/``browser``/``#`` lines are dropped) and converted with
``bedGraphToBigWig``. Use ``--clip`` to clamp values that marginally
overrun the reference bounds. Requires the UCSC kent tools (provided by
the ``lncbookv3-kent`` Docker image).
"""
import argparse
import os
import tempfile

from genomic_intervals.tooling import require_tool, run_command
from genomic_intervals.track import read_chrom_sizes, sort_track_file


def bedgraph_to_bigwig(input_path: str, chrom_sizes: str, output_bigwig: str,
                       clip: bool = False, report_path: str = None):
    order, _sizes = read_chrom_sizes(chrom_sizes)
    converter = require_tool("bedGraphToBigWig")

    fd, sorted_path = tempfile.mkstemp(suffix=".bedgraph")
    os.close(fd)
    try:
        sort_track_file(input_path, sorted_path, order, min_columns=4)
        command = [converter, sorted_path, chrom_sizes, output_bigwig]
        if clip:
            command.append("-clip")
        run_command(command)
    finally:
        os.remove(sorted_path)

    if report_path:
        with open(report_path, "w", encoding="utf-8") as rep:
            rep.write("metric\tvalue\n")
            rep.write("tool\tbedGraphToBigWig\n")
            rep.write(f"clip\t{1 if clip else 0}\n")
            if os.path.exists(output_bigwig):
                rep.write(f"output_bytes\t{os.path.getsize(output_bigwig)}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert bedGraph to bigWig")
    parser.add_argument("--input-bedgraph", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--output-bigwig", required=True)
    parser.add_argument("--clip", action="store_true",
                        help="clip regions outside chromosome bounds")
    parser.add_argument("--report")
    args = parser.parse_args()
    bedgraph_to_bigwig(args.input_bedgraph, args.chrom_sizes, args.output_bigwig,
                       clip=args.clip, report_path=args.report)


if __name__ == "__main__":
    main()
