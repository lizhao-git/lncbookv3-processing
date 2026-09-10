#!/usr/bin/env python3
"""Sort VCF records by contig order and position.

Contig order source (first match wins): ``--contig-order`` (comma-separated
names), the contig names in a ``--fasta-index`` (.fai) file, the ``##contig``
header lines, then first-appearance order for any contigs not yet ranked.
Records on the same contig are ordered by POS while keeping the original
order for ties. All records are held in memory, so this targets typical
site-level VCFs rather than raw whole-genome call sets.
"""
import argparse

from genomic_intervals.io import open_text, open_text_write
from genomic_intervals.vcf import iter_vcf, parse_structured_meta


def load_contig_order(meta_lines, contig_order_text, fasta_index_path):
    if contig_order_text:
        names = [name.strip() for name in contig_order_text.split(",") if name.strip()]
        return names, "explicit"
    if fasta_index_path:
        names = []
        with open_text(fasta_index_path) as (fh, _compression):
            for raw in fh:
                line = raw.strip()
                if line:
                    names.append(line.split("\t")[0])
        return names, "fasta_index"
    names = []
    for line in meta_lines:
        parsed = parse_structured_meta(line)
        if parsed and parsed[0] == "contig":
            name = parsed[1].get("ID")
            if name:
                names.append(name)
    return names, "header"


def sort_vcf(input_path: str, output_vcf: str, report_path=None,
             contig_order_text=None, fasta_index_path=None):
    meta_lines = []
    header_line = None
    records = []
    for _line_no, kind, item in iter_vcf(input_path):
        if kind == "meta":
            meta_lines.append(item)
        elif kind == "header":
            header_line = item
        else:
            records.append(item)

    names, order_source = load_contig_order(meta_lines, contig_order_text, fasta_index_path)
    rank = {}
    for index, name in enumerate(names):
        rank.setdefault(name, index)
    next_rank = len(rank)
    for record in records:
        if record.chrom not in rank:
            rank[record.chrom] = next_rank
            next_rank += 1

    records.sort(key=lambda record: (rank[record.chrom], record.pos))

    with open_text_write(output_vcf) as out:
        for line in meta_lines:
            out.write(line + "\n")
        if header_line:
            out.write(header_line + "\n")
        for record in records:
            out.write(record.to_line() + "\n")

    stats = {
        "input_records": len(records),
        "contigs": len(rank),
        "contig_order_source": order_source,
    }
    if report_path:
        with open(report_path, "w", encoding="utf-8") as rep:
            rep.write("metric\tvalue\n")
            for key, value in stats.items():
                rep.write(f"{key}\t{value}\n")
    return stats


def main():
    parser = argparse.ArgumentParser(description="Sort VCF records by contig order and position")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--contig-order", help="comma-separated contig names in the desired order")
    parser.add_argument("--fasta-index", help=".fai index whose contig order should be used")
    parser.add_argument("--report")
    args = parser.parse_args()
    sort_vcf(args.input_vcf, args.output_vcf, report_path=args.report,
             contig_order_text=args.contig_order, fasta_index_path=args.fasta_index)


if __name__ == "__main__":
    main()
