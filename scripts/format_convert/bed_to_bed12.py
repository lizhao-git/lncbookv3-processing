#!/usr/bin/env python3
"""BED -> BED12 converter.

Converts a BED3-BED6 file to BED12 format by treating each interval as a
single-exon feature. This is useful when a standard BED needs to be passed
to tools that expect BED12 input (e.g. bedToGenePred, bedToBigBed with
BED12 type).

Usage:
    bed_to_bed12.py --input-bed INPUT.bed --output-bed OUTPUT.bed12.bed
                    [--rgb 0]
"""

import argparse
import sys
from pathlib import Path

DESCRIPTION = """\
BED -> BED12 converter.

Each BED3-BED6 interval is expanded to a BED12 single-exon record.
The original 3-6 columns are preserved; the remaining 6 BED12 columns
are filled as:

  thickStart  = start
  thickEnd    = end
  itemRgb     = caller-specified (default "0")
  blockCount  = 1
  blockSizes  = <interval length>,
  blockStarts = 0,
"""


def _die(msg, code=1):
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(code)


def convert_bed_to_bed12(src, dst, rgb):
    """Convert BED to BED12 (single-exon), returning per-record metrics."""
    metrics = {
        "records": 0,
        "strands": {"+": 0, "-": 0, ".": 0},
        "skipped_lines": 0,
    }

    with src.open("rt") as fin, dst.open("wt") as fout:
        for lineno, line in enumerate(fin, 1):
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                metrics["skipped_lines"] += 1
                continue
            fields = stripped.split("\t")
            if len(fields) < 3:
                metrics["skipped_lines"] += 1
                continue

            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                metrics["skipped_lines"] += 1
                continue

            name = fields[3] if len(fields) > 3 else "."
            score = fields[4] if len(fields) > 4 else "0"
            strand = fields[5] if len(fields) > 5 else "."

            block_size = end - start
            block_sizes = f"{block_size},"
            block_starts = "0,"

            fout.write(
                f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t"
                f"{start}\t{end}\t{rgb}\t1\t{block_sizes}\t{block_starts}\n"
            )

            metrics["records"] += 1
            if strand in metrics["strands"]:
                metrics["strands"][strand] += 1
            else:
                metrics["strands"]["."] += 1

    return metrics


def main():
    ap = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--input-bed", required=True, type=Path,
                    help="Input BED file (BED3-BED6+).")
    ap.add_argument("--output-bed", required=True, type=Path,
                    help="Output BED12 file path.")
    ap.add_argument("--rgb", default="0",
                    help="itemRgb value (default: 0)")
    args = ap.parse_args()

    if not args.input_bed.is_file():
        _die(f"input not found: {args.input_bed}")

    args.output_bed.parent.mkdir(parents=True, exist_ok=True)

    metrics = convert_bed_to_bed12(args.input_bed, args.output_bed, args.rgb)

    print(f"records\t{metrics['records']}", file=sys.stderr)
    print(f"skipped_lines\t{metrics['skipped_lines']}", file=sys.stderr)
    for strand, count in sorted(metrics["strands"].items()):
        if count:
            print(f"strand:{strand}\t{count}", file=sys.stderr)

    if metrics["records"] == 0:
        _die("no valid BED records found in input")


if __name__ == "__main__":
    main()
