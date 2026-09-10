#!/usr/bin/env python3
import argparse
import gzip
import io
import sys
import zipfile
from contextlib import contextmanager


def detect_compression(path: str) -> str:
    if path.lower().endswith(".gz"):
        return "gzip"
    if path.lower().endswith(".zip"):
        return "zip"
    try:
        with open(path, "rb") as fh:
            signature = fh.read(4)
            if signature[:2] == b"\x1f\x8b":
                return "gzip"
            if signature == b"PK\x03\x04":
                return "zip"
    except OSError:
        return "plain"
    return "plain"


@contextmanager
def open_text(path: str):
    compression = detect_compression(path)
    if compression == "gzip":
        fh = gzip.open(path, "rt", encoding="utf-8", errors="replace")
        try:
            yield fh, compression
        finally:
            fh.close()
        return

    if compression == "zip":
        zf = zipfile.ZipFile(path, "r")
        members = [name for name in zf.namelist() if not name.endswith("/")]
        if not members:
            zf.close()
            raise SystemExit(f"ZIP archive is empty: {path}")
        vcf_members = [name for name in members if name.lower().endswith(".vcf")]
        member = vcf_members[0] if vcf_members else members[0]
        raw = zf.open(member, "r")
        fh = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
        try:
            yield fh, compression
        finally:
            fh.close()
            raw.close()
            zf.close()
        return

    fh = open(path, "r", encoding="utf-8", errors="replace")
    try:
        yield fh, compression
    finally:
        fh.close()


def main():
    parser = argparse.ArgumentParser(description="Basic ClinVar VCF format validation")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    has_fileformat = False
    has_header = False
    has_clnsig_header = False
    record_count = 0
    errors = []
    warnings = []
    compression = "plain"

    with open_text(args.input_vcf) as (fh, compression):
        for ln, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue

            if line.startswith("##"):
                if line.startswith("##fileformat=VCF"):
                    has_fileformat = True
                if line.startswith("##INFO=<ID=CLNSIG"):
                    has_clnsig_header = True
                continue

            if line.startswith("#CHROM"):
                cols = line.split("\t")
                if len(cols) < 8:
                    errors.append((ln, "#CHROM header has fewer than 8 columns"))
                has_header = True
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                errors.append((ln, f"Expected >= 8 columns, got {len(cols)}"))
                continue

            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
            try:
                pos_i = int(pos)
                if pos_i <= 0:
                    errors.append((ln, "POS must be > 0"))
            except ValueError:
                errors.append((ln, "POS is not an integer"))

            if ref in {"", "."}:
                errors.append((ln, "REF is empty or '.'"))
            if alt in {"", "."}:
                warnings.append((ln, "ALT is empty or '.' (allowed for some ClinVar records)"))

            record_count += 1

    if not has_fileformat:
        errors.append((0, "Missing ##fileformat=VCF header"))
    if not has_header:
        errors.append((0, "Missing #CHROM header line"))
    if not has_clnsig_header:
        errors.append((0, "Missing ##INFO=<ID=CLNSIG...> header (recommended for ClinVar VCF)"))

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"records\t{record_count}\n")
        rep.write(f"has_fileformat_header\t{has_fileformat}\n")
        rep.write(f"has_chrom_header\t{has_header}\n")
        rep.write(f"has_clnsig_header\t{has_clnsig_header}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join([f"line {ln}: {msg}" for ln, msg in errors[:100]]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join([f"line {ln}: {msg}" for ln, msg in warnings[:100]]) + "\n")

    if errors:
        sys.stderr.write("VCF validation failed. See report for details.\n")
        sys.exit(1)

    with open_text(args.input_vcf) as (in_fh, _), open(args.output_vcf, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)


if __name__ == "__main__":
    main()
