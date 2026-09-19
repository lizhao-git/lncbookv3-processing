#!/usr/bin/env python3
"""Normalize an arbitrary BED file into the 10-column feature-interval schema
used by the annotation pipelines (same schema as ``extract_annotation_features``).

Output columns (0-based BED coordinates are kept as-is)::

    chrom  start  end  feature_type  feature_id  gene_id  transcript_id  strand  gene_type  transcript_type

Mapping rules for the optional BED ``name`` column (column 4):

* recognised feature names (``exon``, ``CDS``, ``5UTR``, ``five_prime_utr``, ...)
  are canonicalised into the standard feature-type vocabulary;
* any other name becomes the ``feature_id`` while the row is typed as ``interval``;
* a missing name yields ``feature_type=interval`` and a coordinate-derived id.

``gene_id`` / ``transcript_id`` / ``gene_type`` / ``transcript_type`` cannot be
derived from a plain BED, so they are filled with ``NA``. The strand is taken
from column 6 when it is a valid ``+``/``-`` value, otherwise ``.``.
"""
import argparse
import os
from collections import Counter

from genomic_intervals.io import open_text

# canonicalisation table for common feature-name spellings
CANONICAL_TYPES = {
    "gene": "gene",
    "transcript": "transcript",
    "mrna": "transcript",
    "lnc_rna": "transcript",
    "lncrna": "transcript",
    "exon": "exon",
    "intron": "intron",
    "cds": "cds",
    "utr": "utr",
    "utr5": "utr5",
    "5utr": "utr5",
    "five_prime_utr": "utr5",
    "5_prime_utr": "utr5",
    "utr3": "utr3",
    "3utr": "utr3",
    "three_prime_utr": "utr3",
    "3_prime_utr": "utr3",
    "start_codon": "start_codon",
    "stop_codon": "stop_codon",
}


def normalize_bed(input_path: str, output_path: str, report_path: str, feature_type: str | None = None):
    errors = []
    type_counter = Counter()
    records = 0

    # Rows are streamed to a temporary file while scanning so memory stays
    # constant for large BED inputs; the temp file is renamed to the output
    # path only after validation succeeds.
    tmp_path = f"{output_path}.tmp"
    with open_text(input_path, preferred_exts=(".bed",)) as (fh, _compression), open(
        tmp_path, "w", encoding="utf-8"
    ) as out_fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            records += 1
            cols = line.split()
            if len(cols) < 3:
                errors.append((line_no, f"Expected >= 3 columns, got {len(cols)}"))
                continue
            chrom, start, end = cols[0], cols[1], cols[2]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                errors.append((line_no, "Start/End are not integers"))
                continue
            if start_i < 0:
                errors.append((line_no, "BED start must be >= 0"))
            if end_i <= start_i:
                errors.append((line_no, "BED end must be > start"))
                continue

            name = cols[3] if len(cols) >= 4 else ""
            strand = cols[5] if len(cols) >= 6 and cols[5] in {"+", "-"} else "."

            if feature_type:
                ftype = feature_type
            else:
                ftype = CANONICAL_TYPES.get(name.lower(), name if name else "interval")
            if name and name.lower() not in CANONICAL_TYPES:
                fid = name
            else:
                fid = f"{chrom}:{start_i + 1}-{end_i}"

            type_counter[ftype] += 1
            out_fh.write(
                "\t".join([chrom, str(start_i), str(end_i), ftype, fid, "NA", "NA", strand, "NA", "NA"]) + "\n"
            )

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"records\t{records}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        for ftype, count in sorted(type_counter.items()):
            rep.write(f"type:{ftype}\t{count}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")

    if errors:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise SystemExit("BED normalization failed. See report for details.")

    os.replace(tmp_path, output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Normalize a BED file into the 10-column feature-interval schema"
    )
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument(
        "--feature-type",
        help="Force every row to this feature type instead of deriving it from the BED name column",
    )
    args = parser.parse_args()
    normalize_bed(args.input_bed, args.output_bed, args.report, feature_type=args.feature_type)


if __name__ == "__main__":
    main()
