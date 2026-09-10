#!/usr/bin/env python3
"""Split a VCF into per-contig or per-sample files.

``--mode chrom`` writes one VCF per contig and ``--mode sample`` one VCF
per sample column, all into ``--out-dir`` named
``PREFIX.<key>.vcf[.gz]`` (keys are sanitized for file names). The output
compression mirrors the input. Every output keeps the full meta header;
per-contig outputs keep all samples, per-sample outputs keep only the
selected sample column.
"""
import argparse
import os
import re

from genomic_intervals.io import open_text_write
from genomic_intervals.vcf import header_samples, iter_vcf, build_header_line

SANITIZE_RE = re.compile(r"[^A-Za-z0-9._-]")


def sanitize(name: str) -> str:
    return SANITIZE_RE.sub("_", name)


def split_by_chrom(input_path: str, out_dir: str, prefix: str, suffix: str):
    meta_lines = []
    header_line = None
    handles = {}
    file_names = {}
    counts = {}
    try:
        for _line_no, kind, item in iter_vcf(input_path):
            if kind == "meta":
                meta_lines.append(item)
                continue
            if kind == "header":
                header_line = item
                continue
            chrom = item.chrom
            if chrom not in handles:
                file_name = f"{prefix}.{sanitize(chrom)}{suffix}"
                file_names[chrom] = file_name
                handle = open_text_write(os.path.join(out_dir, file_name))
                fh = handle.__enter__()
                for line in meta_lines:
                    fh.write(line + "\n")
                if header_line:
                    fh.write(header_line + "\n")
                handles[chrom] = (handle, fh)
                counts[chrom] = 0
            handles[chrom][1].write(item.to_line() + "\n")
            counts[chrom] += 1
    finally:
        for handle, _fh in handles.values():
            handle.__exit__(None, None, None)
    return {file_names[key]: counts[key] for key in counts}


def split_by_sample(input_path: str, out_dir: str, prefix: str, suffix: str):
    meta_lines = []
    samples = []
    handles = {}
    counts = {}
    try:
        for _line_no, kind, item in iter_vcf(input_path):
            if kind == "meta":
                meta_lines.append(item)
                continue
            if kind == "header":
                samples = header_samples(item)
                if not samples:
                    raise SystemExit("VCF contains no sample columns to split on.")
                for sample in samples:
                    file_name = f"{prefix}.{sanitize(sample)}{suffix}"
                    handle = open_text_write(os.path.join(out_dir, file_name))
                    fh = handle.__enter__()
                    for line in meta_lines:
                        fh.write(line + "\n")
                    fh.write(build_header_line([sample]) + "\n")
                    handles[sample] = (handle, fh, file_name)
                    counts[sample] = 0
                continue
            if len(item.samples) != len(samples):
                raise SystemExit(
                    f"Record at line {item.line_no} has {len(item.samples)} sample columns "
                    f"but the header declares {len(samples)}."
                )
            if not item.formats:
                raise SystemExit(f"Record at line {item.line_no} has no FORMAT column.")
            for index, sample in enumerate(samples):
                cols = [item.chrom, str(item.pos), item.id, item.ref, item.alt_text(),
                        item.qual, item.filt, item.info, ":".join(item.formats),
                        item.samples[index]]
                handles[sample][1].write("\t".join(cols) + "\n")
                counts[sample] += 1
    finally:
        for handle, _fh, _name in handles.values():
            handle.__exit__(None, None, None)
    return {handles[sample][2]: counts[sample] for sample in counts}


def main():
    parser = argparse.ArgumentParser(description="Split a VCF by chromosome or by sample")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--mode", choices=["chrom", "sample"], required=True)
    parser.add_argument("--prefix", default="records")
    parser.add_argument("--report")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    suffix = ".vcf.gz" if args.input_vcf.lower().endswith(".gz") else ".vcf"
    if args.mode == "chrom":
        results = split_by_chrom(args.input_vcf, args.out_dir, args.prefix, suffix)
    else:
        results = split_by_sample(args.input_vcf, args.out_dir, args.prefix, suffix)

    if args.report:
        with open(args.report, "w", encoding="utf-8") as rep:
            rep.write("metric\tvalue\n")
            rep.write(f"mode\t{args.mode}\n")
            rep.write(f"output_files\t{len(results)}\n")
            for name in sorted(results):
                rep.write(f"records_{name}\t{results[name]}\n")


if __name__ == "__main__":
    main()
