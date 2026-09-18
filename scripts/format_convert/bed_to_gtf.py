#!/usr/bin/env python3
"""BED -> GTF converter.

Converts a BED file (BED3-BED6+) to GTF format. Each BED interval becomes
an 'exon' feature in the GTF output. The BED name column (field 4) is used
as both gene_id and transcript_id in the GTF attributes. The start
coordinate is converted from 0-based (BED) to 1-based (GTF).

Usage:
    bed_to_gtf.py --input-bed INPUT.bed --output-gtf OUTPUT.gtf
                  [--source bed_to_gtf] [--feature exon]
"""

import argparse
import sys
from pathlib import Path

DESCRIPTION = """\
BED -> GTF converter.

Each BED interval becomes a GTF feature line. The BED name column
(field 4) is used as both gene_id and transcript_id in the attributes.
The start coordinate is converted from 0-based (BED) to 1-based (GTF).

Output GTF columns:
  seqname     = BED chrom
  source      = caller-specified (default "bed_to_gtf")
  feature     = caller-specified (default "exon")
  start       = BED start + 1  (1-based, GTF convention)
  end         = BED end
  score       = BED score (field 5, or ".")
  strand      = BED strand (field 6, or ".")
  frame       = "."
  attributes  = gene_id "NAME"; transcript_id "NAME";
"""


def _die(msg, code=1):
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(code)


def convert_bed_to_gtf(src, dst, source, feature):
    """Convert BED to GTF, returning per-record metrics."""
    metrics = {
        "records": 0,
        "strands": {"+": 0, "-": 0, ".": 0},
        "skipped_lines": 0,
        "chroms": set(),
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

            name = fields[3] if len(fields) > 3 else f"feature_{metrics['records'] + 1}"
            score = fields[4] if len(fields) > 4 else "."
            strand = fields[5] if len(fields) > 5 else "."

            # GTF is 1-based
            gtf_start = start + 1
            gtf_end = end

            attrs = f'gene_id "{name}"; transcript_id "{name}";'
            fout.write(
                f"{chrom}\t{source}\t{feature}\t"
                f"{gtf_start}\t{gtf_end}\t{score}\t{strand}\t.\t{attrs}\n"
            )

            metrics["records"] += 1
            metrics["chroms"].add(chrom)
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
    ap.add_argument("--output-gtf", required=True, type=Path,
                    help="Output GTF file path.")
    ap.add_argument("--source", default="bed_to_gtf",
                    help="GTF source column value (default: bed_to_gtf)")
    ap.add_argument("--feature", default="exon",
                    help="GTF feature type (default: exon)")
    args = ap.parse_args()

    if not args.input_bed.is_file():
        _die(f"input not found: {args.input_bed}")

    args.output_gtf.parent.mkdir(parents=True, exist_ok=True)

    metrics = convert_bed_to_gtf(
        args.input_bed, args.output_gtf, args.source, args.feature
    )

    print(f"records\t{metrics['records']}", file=sys.stderr)
    print(f"skipped_lines\t{metrics['skipped_lines']}", file=sys.stderr)
    print(f"chroms\t{len(metrics['chroms'])}", file=sys.stderr)
    for strand, count in sorted(metrics["strands"].items()):
        if count:
            print(f"strand:{strand}\t{count}", file=sys.stderr)

    if metrics["records"] == 0:
        _die("no valid BED records found in input")


if __name__ == "__main__":
    main()
