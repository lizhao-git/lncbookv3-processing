"""Validate FASTA sequence files (stdlib only)."""
from .io import open_text

IUPAC_DNA = set("ACGTURYSWKMBDHVN")
ALLOWED_EXTRA = set("-.*")


def _sequence_state():
    return {"headers": 0, "bases": 0, "min_length": None, "max_length": 0, "empty": 0}


def validate_fasta(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    compression = "plain"
    seen = set()
    duplicates = 0
    invalid_chars = {}
    invalid_lines = 0
    current_header = None
    current_length = 0
    stats = _sequence_state()

    def close_sequence():
        nonlocal current_length
        if current_header is None:
            return
        stats["bases"] += current_length
        if current_length == 0:
            stats["empty"] += 1
        else:
            if stats["min_length"] is None or current_length < stats["min_length"]:
                stats["min_length"] = current_length
            if current_length > stats["max_length"]:
                stats["max_length"] = current_length
        current_length = 0

    with open_text(input_path, preferred_exts=(".fa", ".fasta", ".fna", ".fas")) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                close_sequence()
                header = line[1:].strip()
                if not header:
                    errors.append((line_no, "FASTA header is empty"))
                    current_header = f"line{line_no}"
                else:
                    current_header = header
                stats["headers"] += 1
                if current_header in seen:
                    duplicates += 1
                    warnings.append((line_no, f"Duplicate sequence header: {current_header}"))
                seen.add(current_header)
                continue
            if current_header is None:
                errors.append((line_no, "Sequence data appears before any FASTA header"))
                current_header = f"line{line_no}"
            sequence = line.upper()
            invalid = False
            for char in sequence:
                if char in IUPAC_DNA or char in ALLOWED_EXTRA:
                    continue
                invalid_chars[char] = invalid_chars.get(char, 0) + 1
                invalid = True
            if invalid:
                invalid_lines += 1
                if invalid_lines <= 20:
                    warnings.append((line_no, f"Non-IUPAC character(s) in sequence line: {line[:60]}"))
            current_length += len(sequence)

    close_sequence()

    if stats["headers"] == 0:
        errors.append((0, "No FASTA headers found"))
    if invalid_chars:
        summary = ", ".join(f"{char!r} x{count}" for char, count in sorted(invalid_chars.items()))
        warnings.append((0, f"Non-IUPAC characters present: {summary}"))
    if stats["empty"]:
        warnings.append((0, f"{stats['empty']} sequence(s) have no residues"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"sequences\t{stats['headers']}\n")
        rep.write(f"total_bases\t{stats['bases']}\n")
        rep.write(f"min_length\t{stats['min_length'] if stats['min_length'] is not None else 0}\n")
        rep.write(f"max_length\t{stats['max_length']}\n")
        rep.write(f"empty_sequences\t{stats['empty']}\n")
        rep.write(f"duplicate_headers\t{duplicates}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("FASTA validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".fa", ".fasta", ".fna", ".fas")) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
