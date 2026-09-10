#!/usr/bin/env python3
"""Convert SAM/BAM alignments to BED intervals.

Works for text SAM/SAM.gz directly and pipes BAM through ``samtools view``
when the input name ends with ``.bam``. Each alignment becomes a BED6 row
covering its reference span (CIGAR operations M/D/N/=/X). With
``--split-blocks`` instead of one span per alignment, every aligned block
(M/=/X runs separated by N) is written as its own BED row, which is the
usual choice for spliced RNA alignments. Unmapped reads are skipped.
"""
import argparse
import re
import subprocess

from genomic_intervals.io import open_text
from genomic_intervals.tooling import require_tool

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")
BLOCK_OPS = set("M=X")


def iter_alignment_lines(path: str):
    """Yield SAM text lines from a SAM/SAM.gz file or a samtools view pipe."""
    if path.lower().endswith(".bam"):
        samtools = require_tool("samtools")
        proc = subprocess.Popen(
            [samtools, "view", "--no-PG", path],
            stdout=subprocess.PIPE, text=True,
        )
        try:
            for raw in proc.stdout:
                yield raw.rstrip("\n")
        finally:
            proc.stdout.close()
            if proc.wait() != 0:
                raise SystemExit(f"samtools view failed for {path}")
    else:
        with open_text(path) as (fh, _compression):
            for raw in fh:
                yield raw.rstrip("\n")


def parse_cigar(cigar: str):
    """Return [(length, op)] for a complete CIGAR string, else None."""
    if cigar in ("", "*"):
        return None
    ops = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if not ops or "".join(f"{length}{op}" for length, op in ops) != cigar:
        return None
    return ops


def alignment_to_bed(input_path: str, output_bed: str, min_mapq=None,
                     split_blocks: bool = False, report_path=None):
    stats = {
        "input_lines": 0,
        "header_lines": 0,
        "unmapped_skipped": 0,
        "low_mapq_skipped": 0,
        "malformed_skipped": 0,
        "records_emitted": 0,
    }
    with open(output_bed, "w", encoding="utf-8") as out:
        for line in iter_alignment_lines(input_path):
            if line.startswith("@"):
                stats["header_lines"] += 1
                continue
            stats["input_lines"] += 1
            cols = line.split("\t")
            if len(cols) < 11:
                stats["malformed_skipped"] += 1
                continue
            qname, flag_text, chrom, pos_text, mapq_text, cigar = cols[:6]
            try:
                flag, pos, mapq = int(flag_text), int(pos_text), int(mapq_text)
            except ValueError:
                stats["malformed_skipped"] += 1
                continue
            if flag & 0x4 or chrom == "*" or pos == 0:
                stats["unmapped_skipped"] += 1
                continue
            if min_mapq is not None and mapq < min_mapq:
                stats["low_mapq_skipped"] += 1
                continue
            ops = parse_cigar(cigar)
            if ops is None:
                stats["malformed_skipped"] += 1
                continue

            strand = "-" if flag & 0x10 else "+"
            start = pos - 1
            if split_blocks:
                rows = []
                offset = 0
                for length, op in ops:
                    if op in BLOCK_OPS and length > 0:
                        rows.append((start + offset, start + offset + length))
                    if op in REF_CONSUMING:
                        offset += length
                if not rows:
                    stats["malformed_skipped"] += 1
                    continue
            else:
                span = sum(length for length, op in ops if op in REF_CONSUMING)
                if span <= 0:
                    stats["malformed_skipped"] += 1
                    continue
                rows = [(start, start + span)]

            for row_start, row_end in rows:
                out.write("\t".join([
                    chrom, str(row_start), str(row_end), qname, str(mapq), strand,
                ]) + "\n")
                stats["records_emitted"] += 1

    if report_path:
        with open(report_path, "w", encoding="utf-8") as rep:
            rep.write("metric\tvalue\n")
            for key, value in stats.items():
                rep.write(f"{key}\t{value}\n")
    return stats


def main():
    parser = argparse.ArgumentParser(description="Convert SAM/BAM alignments to BED intervals")
    parser.add_argument("--input-alignment", required=True, help="SAM, SAM.gz or BAM file")
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--min-mapq", type=int, default=None)
    parser.add_argument("--split-blocks", action="store_true",
                        help="emit one BED row per aligned block (M/=/X runs)")
    parser.add_argument("--report")
    args = parser.parse_args()
    alignment_to_bed(args.input_alignment, args.output_bed, min_mapq=args.min_mapq,
                     split_blocks=args.split_blocks, report_path=args.report)


if __name__ == "__main__":
    main()
