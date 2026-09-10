#!/usr/bin/env python3
"""Normalize a VCF by splitting multiallelic records into biallelic ones.

Biallelic records pass through unchanged. For a multiallelic record, one
output line is emitted per ALT allele with the ALT column reduced to that
allele, the GT recoded (kept ALT -> 1, other ALTs -> missing), INFO fields
subset by their declared Number (A: per-ALT element, R: REF + kept ALT),
and FORMAT fields subset the same way plus diploid Number=G subsetting.
Fields whose Number cannot be resolved are kept unchanged and counted in
the report. Left-alignment/trimming requires a reference FASTA and is out
of scope here.
"""
import argparse

from genomic_intervals.io import open_text_write
from genomic_intervals.vcf import (
    format_info,
    iter_vcf,
    parse_info,
    parse_structured_meta,
    recode_gt,
    select_number_value,
    subset_genotype_field,
)


def _split_info(info_text, info_defs, alt_index, alt_count):
    items = []
    for key, value in parse_info(info_text):
        if value is None:
            items.append((key, None))
            continue
        number = info_defs.get(key, {}).get("Number")
        selected = select_number_value(value, number, alt_index, alt_count)
        items.append((key, ",".join(selected)))
    return format_info(items)


def _split_sample(sample_text, formats, format_defs, alt_index, alt_count, stats):
    values = sample_text.split(":")
    if len(values) < len(formats):
        values = values + ["."] * (len(formats) - len(values))
    out = []
    for field, value in zip(formats, values):
        if field == "GT":
            out.append(recode_gt(value, alt_index))
            continue
        if value in (".", ""):
            out.append(value)
            continue
        number = format_defs.get(field, {}).get("Number")
        if number == "G":
            subset = subset_genotype_field(value, alt_count, alt_index)
            if subset is None:
                stats["genotype_fields_kept"] += 1
                out.append(value)
            else:
                out.append(",".join(subset))
            continue
        selected = select_number_value(value, number, alt_index, alt_count)
        out.append(",".join(selected))
    return ":".join(out)


def _split_record_line(record, alt_index, alt, info_defs, format_defs, stats):
    alt_count = len(record.alts)
    info_text = _split_info(record.info, info_defs, alt_index, alt_count)
    cols = [record.chrom, str(record.pos), record.id, record.ref, alt,
            record.qual, record.filt, info_text]
    if record.formats:
        cols.append(":".join(record.formats))
        for sample in record.samples:
            cols.append(_split_sample(sample, record.formats, format_defs,
                                      alt_index, alt_count, stats))
    elif record.samples:
        cols.extend(record.samples)
    return "\t".join(cols)


def normalize_vcf(input_path: str, output_vcf: str, report_path: str = None):
    stats = {
        "input_records": 0,
        "multiallelic_records": 0,
        "output_records": 0,
        "genotype_fields_kept": 0,
    }
    info_defs = {}
    format_defs = {}
    with open_text_write(output_vcf) as out:
        for _line_no, kind, item in iter_vcf(input_path):
            if kind == "meta":
                parsed = parse_structured_meta(item)
                if parsed:
                    meta_kind, attrs = parsed
                    name = attrs.get("ID")
                    if name:
                        if meta_kind == "INFO":
                            info_defs[name] = attrs
                        elif meta_kind == "FORMAT":
                            format_defs[name] = attrs
                out.write(item + "\n")
                continue
            if kind == "header":
                out.write(item + "\n")
                continue

            stats["input_records"] += 1
            if len(item.alts) <= 1:
                out.write(item.to_line() + "\n")
                stats["output_records"] += 1
                continue
            stats["multiallelic_records"] += 1
            for alt_index, alt in enumerate(item.alts):
                out.write(_split_record_line(item, alt_index, alt,
                                             info_defs, format_defs, stats) + "\n")
                stats["output_records"] += 1

    if report_path:
        with open(report_path, "w", encoding="utf-8") as rep:
            rep.write("metric\tvalue\n")
            for key, value in stats.items():
                rep.write(f"{key}\t{value}\n")
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Split multiallelic VCF records into biallelic records")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--report")
    args = parser.parse_args()
    normalize_vcf(args.input_vcf, args.output_vcf, args.report)


if __name__ == "__main__":
    main()
