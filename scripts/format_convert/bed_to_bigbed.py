#!/usr/bin/env python3
"""Convert a BED file into a bigBed track using the kent tool.

The BED type (``bedN``) is inferred from the (uniform) column count of the
data rows unless ``--bed-type`` is given. Rows are sorted into the order
declared by ``--chrom-sizes`` and converted with ``bedToBigBed``. Requires
the UCSC kent tools, which are installed in the lncbookv3-processing
Docker image.
"""
import argparse
import os
import tempfile

from genomic_intervals.io import open_text
from genomic_intervals.tooling import require_tool, run_command
from genomic_intervals.track import read_chrom_sizes, sort_track_file


def infer_bed_type(input_path: str) -> str:
    first = None
    with open_text(input_path) as (fh, _compression):
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            columns = len(line.split("\t"))
            if first is None:
                first = columns
            elif columns != first:
                raise SystemExit(
                    f"{input_path}: line {line_no}: inconsistent column count "
                    f"({columns} vs {first}); bigBed requires a uniform BED type"
                )
    if first is None:
        raise SystemExit(f"{input_path}: no BED rows found")
    if not 3 <= first <= 12:
        raise SystemExit(f"{input_path}: BED column count must be between 3 and 12, got {first}")
    return f"bed{first}"


def bed_to_bigbed(input_path: str, chrom_sizes: str, output_bigbed: str,
                  bed_type: str = None, report_path: str = None):
    order, _sizes = read_chrom_sizes(chrom_sizes)
    converter = require_tool("bedToBigBed")
    if not bed_type:
        bed_type = infer_bed_type(input_path)

    fd, sorted_path = tempfile.mkstemp(suffix=".bed")
    os.close(fd)
    try:
        sort_track_file(input_path, sorted_path, order, min_columns=3)
        run_command([converter, f"-type={bed_type}", sorted_path, chrom_sizes, output_bigbed])
    finally:
        os.remove(sorted_path)

    if report_path:
        with open(report_path, "w", encoding="utf-8") as rep:
            rep.write("metric\tvalue\n")
            rep.write("tool\tbedToBigBed\n")
            rep.write(f"bed_type\t{bed_type}\n")
            if os.path.exists(output_bigbed):
                rep.write(f"output_bytes\t{os.path.getsize(output_bigbed)}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert BED to bigBed")
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--output-bigbed", required=True)
    parser.add_argument("--bed-type", default=None,
                        help="override the inferred bedN type, e.g. bed6")
    parser.add_argument("--report")
    args = parser.parse_args()
    bed_to_bigbed(args.input_bed, args.chrom_sizes, args.output_bigbed,
                  bed_type=args.bed_type, report_path=args.report)


if __name__ == "__main__":
    main()
