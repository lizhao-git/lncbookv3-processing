"""Validate bedMethyl files (stdlib only).

bedMethyl is the tabular base-modification format produced by ``modkit
pileup`` (18 columns, https://software-docs.nanoporetech.com/modkit/intro_pileup.html)
and described by ENCODE (11 columns, no per-class counts). Both layouts are
accepted; arithmetic checks are applied to the full modkit layout:

- ``Nvalid_cov = Nmod + Nother_mod + Ncanonical``
- ``percent_modified = Nmod / Nvalid_cov * 100``
"""
import re

from .io import open_text

MODKIT_COLUMNS = 18
ENCODE_COLUMNS = 11

# Column 4 is a single-letter mod code, a ChEBI number, or "code,motif,offset"
# when modkit pileup is run with more than one motif.
MOD_CODE_RE = re.compile(r"^([A-Za-z]+|\d+)(,[A-Za-z]+,\d+)?$")
COLOR_RE = re.compile(r"^(\d{1,3}),(\d{1,3}),(\d{1,3})$")

PERCENT_TOLERANCE = 0.01


def _parse_int(errors, line_no, label, token, minimum=0):
    """Parse a non-negative integer field, recording errors and returning None."""
    try:
        value = int(token)
    except ValueError:
        errors.append((line_no, f"{label} is not an integer: {token!r}"))
        return None
    if value < minimum:
        errors.append((line_no, f"{label} must be >= {minimum}, got {value}"))
        return None
    return value


def validate_bedmethyl(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    compression = "plain"
    records = 0
    layout = None
    chroms = []
    mod_types = {}
    strand_counts = {"+": 0, "-": 0, ".": 0}
    valid_cov_total = 0
    min_cov = None
    max_cov = None
    duplicates = 0
    unsorted = 0
    multi_base = 0
    zero_coverage = 0
    score_mismatch = 0
    percent_mismatch = 0
    seen = {}
    last_chrom = None
    last_start = None

    with open_text(input_path, preferred_exts=(".bedmethyl", ".bed.gz", ".bed", ".txt")) as (fh, detected):
        compression = detected
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) not in (MODKIT_COLUMNS, ENCODE_COLUMNS):
                errors.append(
                    (line_no, f"expected {MODKIT_COLUMNS} (modkit) or {ENCODE_COLUMNS} (ENCODE) "
                              f"tab-separated columns, got {len(fields)}")
                )
                continue
            if layout is None:
                layout = len(fields)
            elif len(fields) != layout:
                errors.append(
                    (line_no, f"column count {len(fields)} differs from earlier records ({layout})")
                )
                continue

            records += 1
            chrom, start, end, mod_code, score, strand = fields[:6]

            if not chrom:
                errors.append((line_no, "chromosome name is empty"))
            if chrom not in chroms:
                chroms.append(chrom)

            start_i = _parse_int(errors, line_no, "start", start, minimum=0)
            end_i = _parse_int(errors, line_no, "end", end, minimum=0)
            if start_i is not None and end_i is not None:
                if end_i <= start_i:
                    errors.append((line_no, f"end ({end_i}) must be > start ({start_i})"))
                elif end_i - start_i > 1:
                    multi_base += 1

            if not MOD_CODE_RE.match(mod_code):
                errors.append(
                    (line_no, f"mod code {mod_code!r} is not a single-letter code, a ChEBI "
                              "number or code,motif,offset")
                )
            else:
                code = mod_code.split(",")[0]
                mod_types[code] = mod_types.get(code, 0) + 1

            if strand not in strand_counts:
                errors.append((line_no, f"strand must be '+', '-' or '.', got {strand!r}"))
            else:
                strand_counts[strand] += 1

            score_i = _parse_int(errors, line_no, "score", score, minimum=0)

            # Columns 7-8 repeat start/end for compatibility; column 9 is RGB.
            thick_start = _parse_int(errors, line_no, "thickStart", fields[6], minimum=0)
            thick_end = _parse_int(errors, line_no, "thickEnd", fields[7], minimum=0)
            if (thick_start is not None and start_i is not None and thick_start != start_i) or \
                    (thick_end is not None and end_i is not None and thick_end != end_i):
                warnings.append((line_no, "thickStart/thickEnd differ from start/end"))

            color = COLOR_RE.match(fields[8])
            if not color or any(int(part) > 255 for part in color.groups()):
                errors.append(
                    (line_no, f"color {fields[8]!r} is not a R,G,B triple with components 0-255")
                )

            cov_i = None
            if layout == MODKIT_COLUMNS:
                counts = {}
                for label, token, idx in (
                    ("Nvalid_cov", fields[9], 9),
                    ("Nmod", fields[11], 11),
                    ("Ncanonical", fields[12], 12),
                    ("Nother_mod", fields[13], 13),
                    ("Ndelete", fields[14], 14),
                    ("Nfail", fields[15], 15),
                    ("Ndiff", fields[16], 16),
                    ("Nnocall", fields[17], 17),
                ):
                    counts[label] = _parse_int(errors, line_no, label, fields[idx], minimum=0)

                cov_i = counts["Nvalid_cov"]
                if cov_i is not None:
                    valid_cov_total += cov_i
                    min_cov = cov_i if min_cov is None else min(min_cov, cov_i)
                    max_cov = cov_i if max_cov is None else max(max_cov, cov_i)
                    if cov_i == 0:
                        zero_coverage += 1
                    if score_i is not None and score_i != cov_i:
                        score_mismatch += 1
                        warnings.append(
                            (line_no, f"score ({score_i}) differs from Nvalid_cov ({cov_i})")
                        )

                components = (counts["Nmod"], counts["Nother_mod"], counts["Ncanonical"])
                if cov_i is not None and all(part is not None for part in components):
                    total = sum(components)
                    if total != cov_i:
                        errors.append(
                            (line_no, f"Nvalid_cov ({cov_i}) != Nmod + Nother_mod + Ncanonical ({total})")
                        )

                pct_token = fields[10]
                try:
                    pct = float(pct_token)
                except ValueError:
                    errors.append((line_no, f"percent modified is not numeric: {pct_token!r}"))
                    pct = None
                if pct is not None:
                    if pct != pct or pct in (float("inf"), float("-inf")):
                        errors.append((line_no, "percent modified is not finite"))
                    elif pct < 0 or pct > 100:
                        errors.append((line_no, f"percent modified must be within 0-100, got {pct}"))
                    elif cov_i and cov_i > 0 and counts["Nmod"] is not None:
                        expected = counts["Nmod"] / cov_i * 100
                        if abs(pct - expected) > PERCENT_TOLERANCE:
                            percent_mismatch += 1
                            errors.append(
                                (line_no, f"percent modified ({pct}) != Nmod / Nvalid_cov "
                                          f"({expected:.2f})")
                            )
            else:
                try:
                    pct = float(fields[10])
                    if pct != pct or pct < 0 or pct > 100:
                        errors.append((line_no, f"percent modified must be within 0-100, got {pct}"))
                except ValueError:
                    errors.append((line_no, f"percent modified is not numeric: {fields[10]!r}"))
                cov_i = _parse_int(errors, line_no, "coverage", fields[9], minimum=0)

            key = (chrom, start, end, strand, mod_code)
            if key in seen:
                duplicates += 1
                warnings.append((line_no, f"duplicate position {chrom}:{start}-{end} {strand} "
                                          f"{mod_code} (first seen on line {seen[key]})"))
            else:
                seen[key] = line_no

            if start_i is not None:
                if last_start is not None:
                    if chrom == last_chrom and start_i < last_start:
                        unsorted += 1
                    elif chrom != last_chrom and chrom in chroms[:-1]:
                        unsorted += 1
                last_chrom, last_start = chrom, start_i

    if records == 0:
        errors.append((0, "No bedMethyl records found"))
    if unsorted:
        warnings.append((0, f"{unsorted} record(s) break chromosome/position ordering"))
    if multi_base:
        warnings.append((0, f"{multi_base} record(s) span more than one base"))
    if zero_coverage:
        warnings.append((0, f"{zero_coverage} record(s) have zero Nvalid_cov"))

    with open(report_path, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"input_compression\t{compression}\n")
        report.write(f"layout\t{'modkit-18' if layout == MODKIT_COLUMNS else 'encode-11' if layout else 'unknown'}\n")
        report.write(f"records\t{records}\n")
        report.write(f"chrom_count\t{len(chroms)}\n")
        report.write(f"mod_type_count\t{len(mod_types)}\n")
        report.write("mod_types\t" + ",".join(f"{code}:{count}" for code, count in sorted(mod_types.items())) + "\n")
        report.write(f"strand_plus\t{strand_counts['+']}\n")
        report.write(f"strand_minus\t{strand_counts['-']}\n")
        report.write(f"strand_dot\t{strand_counts['.']}\n")
        if layout == MODKIT_COLUMNS:
            report.write(f"total_valid_cov\t{valid_cov_total}\n")
            report.write(f"min_valid_cov\t{min_cov if min_cov is not None else 0}\n")
            report.write(f"max_valid_cov\t{max_cov if max_cov is not None else 0}\n")
        report.write(f"duplicate_positions\t{duplicates}\n")
        report.write(f"unsorted_records\t{unsorted}\n")
        report.write(f"multi_base_records\t{multi_base}\n")
        report.write(f"score_mismatch\t{score_mismatch}\n")
        report.write(f"percent_mismatch\t{percent_mismatch}\n")
        report.write(f"error_count\t{len(errors)}\n")
        report.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            report.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            report.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("bedMethyl validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".bedmethyl", ".bed.gz", ".bed", ".txt")) as (in_fh, _), \
            open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
