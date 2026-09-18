#!/usr/bin/env python3
"""bedMethyl → BED converter.

Converts a bedMethyl file (modkit output, >=10 tab-delimited columns) to a
standard BED file. By default outputs BED6 (chrom, start, end, mod_code,
score, strand). Optional flags append coverage or percent-modification columns.

Usage:
    bedmethyl_to_bed.py --input-bedmethyl INPUT.bedMethyl \
                        --output-bed OUTPUT.bed \
                        [--with-coverage] [--with-percent]
"""

import argparse
import shutil
import sys
from pathlib import Path

DESCRIPTION = """\
bedMethyl → BED converter.

Strips the methylation-specific columns from a modkit bedMethyl file,
producing a standard BED file. The 4th column (mod_code) is preserved as
the BED name field, the 5th (mod_score) as the score, and the 6th (strand)
as the strand.

Optional --with-coverage and --with-percent append the modkit coverage
(column 10) and percent modification (column 11) as extra columns.
"""


def _die(msg, code=1):
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(code)


def convert_bedmethyl(src: Path, dst: Path, with_coverage: bool, with_percent: bool):
    """Convert bedMethyl to BED, yielding per-record metrics."""
    metrics = {
        "records": 0,
        "mod_types": {},
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
            if len(fields) < 6:
                metrics["skipped_lines"] += 1
                continue

            chrom, start, end = fields[0], fields[1], fields[2]
            mod_code = fields[3] if len(fields) > 3 else "."
            score = fields[4] if len(fields) > 4 else "0"
            strand = fields[5] if len(fields) > 5 else "."

            out_cols = [chrom, start, end, mod_code, score, strand]

            if with_coverage and len(fields) > 9:
                out_cols.append(fields[9])
            if with_percent and len(fields) > 10:
                out_cols.append(fields[10])

            fout.write("\t".join(out_cols) + "\n")
            metrics["records"] += 1
            metrics["mod_types"][mod_code] = metrics["mod_types"].get(mod_code, 0) + 1
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
    ap.add_argument("--input-bedmethyl", required=True, type=Path,
                    help="Input bedMethyl file (modkit output).")
    ap.add_argument("--output-bed", required=True, type=Path,
                    help="Output BED file path.")
    ap.add_argument("--with-coverage", action="store_true",
                    help="Append the coverage column (column 10) as an extra field.")
    ap.add_argument("--with-percent", action="store_true",
                    help="Append the percent-modification column (column 11) as an extra field.")
    args = ap.parse_args()

    if not args.input_bedmethyl.is_file():
        _die(f"input not found: {args.input_bedmethyl}")

    args.output_bed.parent.mkdir(parents=True, exist_ok=True)

    metrics = convert_bedmethyl(
        args.input_bedmethyl,
        args.output_bed,
        args.with_coverage,
        args.with_percent,
    )

    # Print metrics to stderr (matching project convention).
    print(f"records\t{metrics['records']}", file=sys.stderr)
    print(f"skipped_lines\t{metrics['skipped_lines']}", file=sys.stderr)
    for mod_code, count in sorted(metrics["mod_types"].items()):
        print(f"mod_type:{mod_code}\t{count}", file=sys.stderr)
    for strand, count in sorted(metrics["strands"].items()):
        if count:
            print(f"strand:{strand}\t{count}", file=sys.stderr)

    if metrics["records"] == 0:
        _die("no valid bedMethyl records found in input")


if __name__ == "__main__":
    main()
