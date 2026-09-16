"""Validate WIG (wiggle) coverage files (stdlib only)."""
from .io import open_text


def _parse_declaration(kind, line, line_no, errors):
    """Parse a fixedStep/variableStep declaration line into a settings dict."""
    settings = {}
    for token in line.split()[1:]:
        key, sep, value = token.partition("=")
        if not sep or not key or not value:
            errors.append((line_no, f"Malformed {kind} declaration token: {token}"))
            continue
        settings[key] = value
    if "chrom" not in settings:
        errors.append((line_no, f"{kind} declaration is missing chrom="))
    for int_key in ("start", "step", "span"):
        if int_key in settings:
            try:
                value_i = int(settings[int_key])
            except ValueError:
                errors.append((line_no, f"{kind} {int_key} is not an integer"))
                continue
            if value_i < 1:
                errors.append((line_no, f"{kind} {int_key} must be >= 1"))
    if kind == "fixedStep" and ("start" not in settings or "step" not in settings):
        errors.append((line_no, "fixedStep declaration requires start= and step="))
    return settings


def validate_wig(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    records = 0
    declarations = 0
    chroms = []
    compression = "plain"
    mode = None
    step = None
    last_pos = None
    out_of_order = 0

    with open_text(input_path, preferred_exts=(".wig",)) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            if line.startswith("fixedStep"):
                declarations += 1
                settings = _parse_declaration("fixedStep", line, line_no, errors)
                mode = "fixed"
                step = settings.get("step")
                last_pos = None
                try:
                    last_pos = int(settings.get("start", "1")) - int(step)
                except ValueError:
                    last_pos = None
                if settings.get("chrom") and settings["chrom"] not in chroms:
                    chroms.append(settings["chrom"])
                continue
            if line.startswith("variableStep"):
                declarations += 1
                settings = _parse_declaration("variableStep", line, line_no, errors)
                mode = "variable"
                last_pos = None
                if settings.get("chrom") and settings["chrom"] not in chroms:
                    chroms.append(settings["chrom"])
                continue
            if mode is None:
                errors.append((line_no, "Data line before any fixedStep/variableStep declaration"))
                continue
            if mode == "fixed":
                try:
                    float(line)
                except ValueError:
                    errors.append((line_no, "fixedStep data line must be a single numeric value"))
                    continue
                records += 1
                if last_pos is not None and step is not None:
                    last_pos += int(step)
            else:
                parts = line.split()
                if len(parts) != 2:
                    errors.append((line_no, "variableStep data line must be 'position value'"))
                    continue
                try:
                    position = int(parts[0])
                    float(parts[1])
                except ValueError:
                    errors.append((line_no, "variableStep data line must be 'position value'"))
                    continue
                if position < 1:
                    errors.append((line_no, "WIG positions are 1-based and must be >= 1"))
                if last_pos is not None and position <= last_pos:
                    out_of_order += 1
                last_pos = position
                records += 1

    if declarations == 0 and records == 0:
        errors.append((0, "No WIG declarations or data found"))
    if out_of_order:
        warnings.append((0, f"{out_of_order} variableStep record(s) are out of order"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"declarations\t{declarations}\n")
        rep.write(f"records\t{records}\n")
        rep.write(f"chrom_count\t{len(chroms)}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("WIG validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".wig",)) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
