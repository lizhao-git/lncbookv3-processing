"""Validate SAM alignment files (stdlib only, no samtools dependency)."""
import re

from .io import open_text

CIGAR_RE = re.compile(r"^(\d+[MIDNSHP=X])+$")
CIGAR_QUERY_CONSUMING = set("MIS=X")
SEQ_RE = re.compile(r"^[A-Za-z=.]+$")
OPTIONAL_TYPES = set("AifZHB")
HEADER_TYPES = ("@HD", "@SQ", "@RG", "@PG", "@CO")


def _validate_cigar(cigar, line_no, errors):
    """Return query-consuming length, or None when the CIGAR is invalid."""
    if cigar == "*":
        return None
    if not CIGAR_RE.match(cigar):
        errors.append((line_no, f"Malformed CIGAR: {cigar}"))
        return None
    query_length = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in CIGAR_QUERY_CONSUMING:
            query_length += int(length)
    return query_length


def validate_sam(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    header_lines = 0
    records = 0
    unmapped = 0
    compression = "plain"
    seen_alignment = False

    with open_text(input_path, preferred_exts=(".sam",)) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("@"):
                if seen_alignment:
                    errors.append((line_no, "Header line appears after alignment records"))
                    continue
                if not line.startswith(HEADER_TYPES):
                    errors.append((line_no, f"Unknown header record: {line.split()[0]}"))
                    continue
                header_lines += 1
                if line.startswith("@HD"):
                    if header_lines != 1:
                        errors.append((line_no, "@HD must be the first header line"))
                    if "VN:" not in line:
                        warnings.append((line_no, "@HD is missing VN:"))
                if line.startswith("@SQ"):
                    fields = dict(token.partition(":")[::2] for token in line.split()[1:] if ":" in token)
                    if "SN" not in fields:
                        errors.append((line_no, "@SQ is missing SN"))
                    if "LN" not in fields:
                        errors.append((line_no, "@SQ is missing LN"))
                    else:
                        try:
                            if int(fields["LN"]) < 1:
                                errors.append((line_no, "@SQ LN must be >= 1"))
                        except ValueError:
                            errors.append((line_no, "@SQ LN is not an integer"))
                continue

            seen_alignment = True
            records += 1
            cols = line.split("\t")
            if len(cols) < 11:
                errors.append((line_no, f"SAM record requires >= 11 fields, got {len(cols)}"))
                continue
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = cols[:11]

            for label, value, low, high in (
                ("FLAG", flag, 0, 65535),
                ("POS", pos, 0, None),
                ("MAPQ", mapq, 0, 255),
                ("PNEXT", pnext, 0, None),
            ):
                try:
                    int_value = int(value)
                except ValueError:
                    errors.append((line_no, f"{label} is not an integer: {value}"))
                    continue
                if int_value < low or (high is not None and int_value > high):
                    errors.append((line_no, f"{label} out of range: {value}"))
            try:
                int(tlen)
            except ValueError:
                errors.append((line_no, f"TLEN is not an integer: {tlen}"))
            try:
                flag_i = int(flag)
                if flag_i & 0x4:
                    unmapped += 1
            except ValueError:
                pass

            if rnext not in ("*", "=") and rnext == rname:
                warnings.append((line_no, "RNEXT repeats RNAME instead of '='"))
            query_length = _validate_cigar(cigar, line_no, errors)
            if seq != "*":
                if not SEQ_RE.match(seq):
                    errors.append((line_no, "SEQ contains characters outside [A-Za-z=.]"))
                if query_length is not None and len(seq) != query_length:
                    errors.append((line_no, f"SEQ length {len(seq)} does not match CIGAR query length {query_length}"))
            if qual != "*":
                if seq != "*" and len(qual) != len(seq):
                    errors.append((line_no, f"QUAL length {len(qual)} does not match SEQ length {len(seq)}"))
                for char in qual:
                    if not 33 <= ord(char) <= 126:
                        errors.append((line_no, "QUAL contains characters outside printable ASCII 33-126"))
                        break
            for extra in cols[11:]:
                parts = extra.split(":")
                if len(parts) < 3 or len(parts[0]) != 2:
                    errors.append((line_no, f"Malformed optional field: {extra}"))
                    continue
                if parts[1] not in OPTIONAL_TYPES:
                    errors.append((line_no, f"Unknown optional field type '{parts[1]}' in: {extra}"))

    if not seen_alignment and not header_lines:
        errors.append((0, "No SAM header or alignment records found"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"header_lines\t{header_lines}\n")
        rep.write(f"records\t{records}\n")
        rep.write(f"unmapped_records\t{unmapped}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("SAM validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".sam",)) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
