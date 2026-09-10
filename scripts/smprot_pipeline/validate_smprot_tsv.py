#!/usr/bin/env python3
import argparse
import gzip
import re
import sys
from typing import Dict, List, Optional, Tuple

REQUIRED_GROUPS = {
    "chrom": {"chrom", "chr", "chromosome"},
    "start": {"start", "chromstart", "tx_start", "pos", "begin"},
    "end": {"end", "chromend", "tx_end", "stop"},
    "protein_id": {"protein_id", "smprot_id", "smprotid", "proteinid", "id", "protein"},
}

OPTIONAL_GROUPS = {
    "evidence": {"evidence", "evidence_type", "support", "support_evidence"},
    "source": {"source", "dataset", "source_db", "species"},
}


def detect_compression(path: str) -> str:
    if path.lower().endswith(".gz"):
        return "gzip"
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
            if magic == b"\x1f\x8b":
                return "gzip"
    except OSError:
        pass
    return "plain"


def open_text(path: str):
    compression = detect_compression(path)
    if compression == "gzip":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace"), compression
    return open(path, "r", encoding="utf-8", errors="replace"), compression


def split_line(line: str, delimiter: str) -> List[str]:
    stripped = line.rstrip("\n\r")
    if delimiter == "tab":
        return stripped.split("\t")
    return re.split(r"\s+", stripped.strip())


def detect_delimiter(line: str) -> str:
    return "tab" if "\t" in line else "whitespace"


def normalize(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_")


def resolve_columns(header: List[str], groups: Dict[str, set]) -> Tuple[Dict[str, int], List[str]]:
    normalized = {normalize(col): idx for idx, col in enumerate(header)}
    resolved = {}
    missing = []
    for canonical, aliases in groups.items():
        idx: Optional[int] = None
        for alias in aliases:
            key = normalize(alias)
            if key in normalized:
                idx = normalized[key]
                break
        if idx is None:
            missing.append(canonical)
        else:
            resolved[canonical] = idx
    return resolved, missing


def main():
    parser = argparse.ArgumentParser(description="Validate SmProt tabular coordinates")
    parser.add_argument("--input-file", action="append", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    errors = []
    warnings = []
    record_count = 0
    kept_count = 0
    file_count = 0

    with open(args.output_tsv, "w", encoding="utf-8") as out:
        out.write("chrom\tstart\tend\tprotein_id\tevidence\tsource\torigin_file\n")

        for input_path in args.input_file:
            file_count += 1

            try:
                fh, compression = open_text(input_path)
            except OSError as exc:
                errors.append((0, f"Cannot open input file {input_path}: {exc}"))
                continue

            with fh:
                first_data_line = None
                first_line_no = 0
                for ln, raw in enumerate(fh, start=1):
                    if raw.strip():
                        first_data_line = raw
                        first_line_no = ln
                        break

                if first_data_line is None:
                    errors.append((0, f"Input file is empty: {input_path}"))
                    continue

                delimiter = detect_delimiter(first_data_line)
                header = split_line(first_data_line, delimiter)

                required_resolved, missing_required = resolve_columns(header, REQUIRED_GROUPS)
                optional_resolved, _ = resolve_columns(header, OPTIONAL_GROUPS)

                for col in missing_required:
                    errors.append((0, f"Missing required column for {col} in {input_path}"))

                if missing_required:
                    continue

                warnings.append((0, f"Detected format for {input_path}: compression={compression}, delimiter={delimiter}"))

                for ln, raw in enumerate(fh, start=first_line_no + 1):
                    if not raw.strip():
                        continue
                    row = split_line(raw, delimiter)
                    record_count += 1

                    if len(row) < len(header):
                        errors.append((ln, f"Row has fewer columns than header in {input_path}"))
                        continue

                    try:
                        start = int(row[required_resolved["start"]].strip())
                        end = int(row[required_resolved["end"]].strip())
                        if start < 0:
                            errors.append((ln, f"start must be >= 0 in {input_path}"))
                            continue
                        if end <= start:
                            errors.append((ln, f"end must be > start in {input_path}"))
                            continue
                    except ValueError:
                        errors.append((ln, f"start/end are not integers in {input_path}"))
                        continue

                    chrom = row[required_resolved["chrom"]].strip()
                    protein_id = row[required_resolved["protein_id"]].strip()
                    evidence = row[optional_resolved["evidence"]].strip() if "evidence" in optional_resolved else "NA"
                    source = row[optional_resolved["source"]].strip() if "source" in optional_resolved else "SmProt"

                    if not protein_id:
                        warnings.append((ln, f"protein_id is empty in {input_path}"))
                        continue
                    if not chrom:
                        warnings.append((ln, f"chrom is empty in {input_path}"))
                        continue

                    out.write(
                        f"{chrom}\t{start}\t{end}\t{protein_id}\t{(evidence or 'NA')}\t{(source or 'SmProt')}\t{input_path}\n"
                    )
                    kept_count += 1

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_file_count\t{file_count}\n")
        rep.write(f"records\t{record_count}\n")
        rep.write(f"kept_records\t{kept_count}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join([f"line {ln}: {msg}" for ln, msg in errors[:100]]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join([f"line {ln}: {msg}" for ln, msg in warnings[:100]]) + "\n")

    if errors:
        sys.stderr.write("SmProt TSV validation failed. See report for details.\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
