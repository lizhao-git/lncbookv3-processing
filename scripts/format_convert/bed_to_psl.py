#!/usr/bin/env python3
"""BED -> PSL converter.

Converts a BED file to PSL format. Each BED interval becomes a PSL
alignment record with 100% match (no mismatches, no gaps). The query
length is set to the interval length; target coordinates use the BED
chrom/start/end directly.

A chrom.sizes file is recommended for accurate tSize values; if not
provided, tSize defaults to the maximum end coordinate seen per chrom.

Usage:
    bed_to_psl.py --input-bed INPUT.bed --output-psl OUTPUT.psl
                  [--target-sizes chrom.sizes]
"""

import argparse
import sys
from pathlib import Path

DESCRIPTION = """\
BED -> PSL converter.

Each BED interval becomes a PSL alignment with 100% identity:
  matches     = interval length
  misMatches  = 0
  repMatches  = 0
  nCount      = 0
  all gaps    = 0
  strand      = BED strand (default "+")
  qName       = BED name
  qSize       = interval length  (simple query model)
  qStart      = 0
  qEnd        = interval length
  tName       = BED chrom
  tSize       = from --target-sizes, or max end per chrom
  tStart      = BED start
  tEnd        = BED end
  blockCount  = 1
  blockSizes  = <interval length>,
  qStarts     = 0,
  tStarts     = <BED start>,

The standard psLayout version 3 header is prepended.
"""

PSL_HEADER_LINES = [
    "psLayout version 3",
    "",
    "\t".join([
        "match", "misMatch", "repMatch", "N's",
        "Q gap count", "Q gap bases", "T gap count", "T gap bases",
        "strand", "Q name", "Q size", "Q start", "Q end",
        "T name", "T size", "T start", "T end",
        "block count", "blockSizes", "qStarts", "tStarts",
    ]),
    "-" * 180,
]


def _die(msg, code=1):
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(code)


def load_sizes(path):
    """Load a chrom.sizes file into a dict."""
    sizes = {}
    with path.open("rt") as fin:
        for line in fin:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 2:
                try:
                    sizes[fields[0]] = int(fields[1])
                except ValueError:
                    continue
    return sizes


def _compute_chrom_max(src):
    """First pass: compute per-chrom max end for tSize fallback."""
    chrom_max = {}
    with src.open("rt") as fin:
        for line in fin:
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 3:
                continue
            try:
                end = int(fields[2])
            except ValueError:
                continue
            chrom = fields[0]
            if chrom not in chrom_max or end > chrom_max[chrom]:
                chrom_max[chrom] = end
    return chrom_max


def convert_bed_to_psl(src, dst, target_sizes=None):
    """Convert BED to PSL, returning per-record metrics."""
    metrics = {
        "records": 0,
        "strands": {"+": 0, "-": 0, ".": 0},
        "skipped_lines": 0,
        "chroms": set(),
        "missing_chroms": 0,
    }

    # First pass: compute per-chrom max end if no sizes file
    chrom_max = {}
    if target_sizes is None:
        chrom_max = _compute_chrom_max(src)

    with src.open("rt") as fin, dst.open("wt") as fout:
        for line in PSL_HEADER_LINES:
            fout.write(line + "\n")

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
            strand = fields[5] if len(fields) > 5 else "+"

            block_size = end - start

            if target_sizes is not None and chrom in target_sizes:
                t_size = target_sizes[chrom]
            elif target_sizes is not None and chrom not in target_sizes:
                t_size = end
                metrics["missing_chroms"] += 1
            else:
                t_size = chrom_max.get(chrom, end)

            psl_fields = [
                str(block_size),  # matches
                "0",              # misMatches
                "0",              # repMatches
                "0",              # nCount
                "0",              # qNumInsert
                "0",              # qBaseInsert
                "0",              # tNumInsert
                "0",              # tBaseInsert
                strand,           # strand
                name,             # qName
                str(block_size),  # qSize
                "0",              # qStart
                str(block_size),  # qEnd
                chrom,            # tName
                str(t_size),      # tSize
                str(start),       # tStart
                str(end),         # tEnd
                "1",              # blockCount
                f"{block_size},", # blockSizes
                "0,",             # qStarts
                f"{start},",      # tStarts
            ]
            fout.write("\t".join(psl_fields) + "\n")

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
                    help="Input BED file.")
    ap.add_argument("--output-psl", required=True, type=Path,
                    help="Output PSL file path.")
    ap.add_argument("--target-sizes", type=Path, default=None,
                    help="chrom.sizes file for accurate tSize values; "
                         "if omitted, tSize defaults to max end per chrom.")
    args = ap.parse_args()

    if not args.input_bed.is_file():
        _die(f"input not found: {args.input_bed}")

    target_sizes = None
    if args.target_sizes is not None:
        if not args.target_sizes.is_file():
            _die(f"target sizes not found: {args.target_sizes}")
        target_sizes = load_sizes(args.target_sizes)

    args.output_psl.parent.mkdir(parents=True, exist_ok=True)

    metrics = convert_bed_to_psl(args.input_bed, args.output_psl, target_sizes)

    print(f"records\t{metrics['records']}", file=sys.stderr)
    print(f"skipped_lines\t{metrics['skipped_lines']}", file=sys.stderr)
    print(f"chroms\t{len(metrics['chroms'])}", file=sys.stderr)
    if metrics["missing_chroms"]:
        print(f"missing_chroms\t{metrics['missing_chroms']}", file=sys.stderr)
    for strand, count in sorted(metrics["strands"].items()):
        if count:
            print(f"strand:{strand}\t{count}", file=sys.stderr)

    if metrics["records"] == 0:
        _die("no valid BED records found in input")


if __name__ == "__main__":
    main()
