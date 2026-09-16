"""Validate bedGraph coverage files (stdlib only)."""
from .io import open_text


def validate_bedgraph(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    records = 0
    chroms = []
    unsorted = 0
    compression = "plain"
    last_chrom = None
    last_start = None

    with open_text(input_path, preferred_exts=(".bedgraph", ".bdg")) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            records += 1
            cols = line.split("\t")
            if len(cols) != 4:
                errors.append((line_no, f"Expected 4 columns, got {len(cols)}"))
                continue
            chrom, start, end, value = cols
            if chrom not in chroms:
                chroms.append(chrom)
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                errors.append((line_no, "Start/End are not integers"))
                continue
            if start_i < 0:
                errors.append((line_no, "bedGraph start must be >= 0"))
            if end_i <= start_i:
                errors.append((line_no, "bedGraph end must be > start"))
            try:
                value_f = float(value)
                if value_f != value_f or value_f in (float("inf"), float("-inf")):
                    errors.append((line_no, "bedGraph value is not finite"))
            except ValueError:
                errors.append((line_no, "bedGraph value is not numeric"))
            if last_start is not None:
                if chrom == last_chrom and start_i < last_start:
                    unsorted += 1
                elif chrom != last_chrom and chrom in chroms[:-1]:
                    unsorted += 1
            last_chrom, last_start = chrom, start_i

    if unsorted:
        warnings.append((0, f"{unsorted} record(s) break chromosome/position ordering"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"records\t{records}\n")
        rep.write(f"chrom_count\t{len(chroms)}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("bedGraph validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".bedgraph", ".bdg")) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
