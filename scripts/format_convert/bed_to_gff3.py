#!/usr/bin/env python3
"""Convert a BED interval file to GFF3.

BED columns 4-6, when present, feed the GFF3 ID, score and strand: the
BED name column becomes the feature ID (falling back to
``chrom:start-end`` in 1-based coordinates), the score column is copied
through, and the strand column is used when it is one of +, - or ..
Coordinates are converted from 0-based half-open (BED) to 1-based
inclusive (GFF3). Extra BED columns (blocks, RGB, ...) are ignored.
"""
import argparse

from genomic_intervals.annotation import gff3_escape
from genomic_intervals.io import open_text

BED_STRANDS = {"+", "-", "."}


def convert_bed_to_gff3(input_path: str, output_path: str, feature: str, source: str) -> int:
    records = 0
    with open_text(input_path) as (fh, _compression), open(output_path, "w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                raise SystemExit(f"{input_path}: line {line_no}: expected at least 3 BED columns, got {len(cols)}")
            chrom, start_text, end_text = cols[0], cols[1], cols[2]
            try:
                start, end = int(start_text), int(end_text)
            except ValueError:
                raise SystemExit(f"{input_path}: line {line_no}: start/end are not integers")
            name = cols[3] if len(cols) > 3 and cols[3] not in ("", ".") else None
            score = cols[4] if len(cols) > 4 and cols[4] not in ("", ".") else "."
            strand = cols[5] if len(cols) > 5 and cols[5] in BED_STRANDS else "."

            start_1based = start + 1
            feature_id = name if name else f"{chrom}:{start_1based}-{end}"
            attrs = f"ID={gff3_escape(feature_id)}"
            out.write("\t".join([
                chrom, source, feature, str(start_1based), str(end), score, strand, ".", attrs,
            ]) + "\n")
            records += 1
    return records


def main():
    parser = argparse.ArgumentParser(description="Convert BED intervals to GFF3")
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--output-gff3", required=True)
    parser.add_argument("--feature", default="region", help="GFF3 feature type to emit (default: region)")
    parser.add_argument("--source", default=".", help="GFF3 source column value (default: .)")
    args = parser.parse_args()
    convert_bed_to_gff3(args.input_bed, args.output_gff3, args.feature, args.source)


if __name__ == "__main__":
    main()
