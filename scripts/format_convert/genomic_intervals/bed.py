from collections import Counter

from .io import open_text


def validate_bed(input_path: str, output_path: str, report_path: str, min_columns: int = 3):
    errors = []
    warnings = []
    chrom_counter = Counter()
    records = 0
    compression = "plain"

    with open_text(input_path, preferred_exts=(".bed",)) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            records += 1
            cols = line.split("\t")
            if len(cols) < min_columns:
                errors.append((line_no, f"Expected >= {min_columns} columns, got {len(cols)}"))
                continue
            chrom, start, end = cols[:3]
            chrom_counter[chrom] += 1
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
            if len(cols) != min_columns and len(cols) < 6:
                warnings.append((line_no, "BED has extra columns but fewer than BED6 columns"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"records\t{records}\n")
        rep.write(f"chrom_count\t{len(chrom_counter)}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("BED validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".bed",)) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
