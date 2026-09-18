"""Validate two-column chrom.sizes files consumed by the UCSC tools (stdlib only).

chrom.sizes is the ``samtools faidx`` column pair (``cut -f1,2``) that
``bedGraphToBigWig``, ``wigToBigWig``, ``bedToBigBed``, ``bedtools sort``,
kent ``liftOver`` and ``pslCheck -targetSizes`` all depend on, so a malformed
file silently corrupts every downstream track.
"""
import re

from .io import open_text

UCSC_NAME_LIMIT = 255

FAI_COLUMNS = 5
CHROM_PREFIXES = ("chr", "chrom")

# Canonical UCSC ordering ranks for the non-numeric chromosome names.
NAMED_RANK = {"x": 23, "y": 24, "m": 25, "mt": 25, "w": 26, "z": 27}


def _split_prefix(name):
    """Return (prefix, stem) splitting a leading chr/chrom prefix when present."""
    lower = name.lower()
    for prefix in CHROM_PREFIXES:
        if lower.startswith(prefix):
            return prefix, name[len(prefix):]
    return "", name


def _order_key(name):
    """Ordering key approximating the UCSC canonical chromosome order."""
    prefix, stem = _split_prefix(name)
    lower = stem.lower()
    if lower in NAMED_RANK:
        return (1, NAMED_RANK[lower], name)
    if stem.isdigit():
        return (1, int(stem), name)
    match = re.search(r"\d+", stem)
    if match:
        return (2, int(match.group()), name)
    return (3, 0, name)


def validate_chrom_sizes(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    compression = "plain"
    records = 0
    duplicates = 0
    space_separated = 0
    long_names = 0
    fai_columns = 0
    sizes = []
    order = []
    prefixes = set()
    seen = {}

    with open_text(input_path, preferred_exts=(".sizes", ".chrom.sizes", ".txt")) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                continue

            tab_fields = line.split("\t")
            if len(tab_fields) != 2:
                tokens = line.split()
                if len(tokens) == 2:
                    space_separated += 1
                    tab_fields = tokens
                elif len(tokens) == FAI_COLUMNS:
                    fai_columns += 1
                    errors.append(
                        (line_no, f"line has {FAI_COLUMNS} columns like a samtools .fai index; "
                                  "chrom.sizes needs only the first two (cut -f1,2)")
                    )
                    continue
                else:
                    errors.append((line_no, f"chrom.sizes needs 2 tab-separated columns, "
                                            f"got {len(tab_fields)}"))
                    continue
            else:
                # Tab separated but a stray space inside the name still breaks the kent tools.
                if " " in tab_fields[0] or " " in tab_fields[1]:
                    errors.append((line_no, "column contains a space; chrom.sizes must not be space padded"))
                    continue

            name, size_token = tab_fields
            if not name:
                errors.append((line_no, "chromosome name is empty"))
                continue

            records += 1
            if len(name) > UCSC_NAME_LIMIT:
                long_names += 1
                warnings.append((line_no, f"chromosome name is {len(name)} characters; the UCSC tools "
                                          f"limit names to {UCSC_NAME_LIMIT}"))
            prefix, _ = _split_prefix(name)
            prefixes.add(prefix)
            order.append(_order_key(name))

            if name in seen:
                duplicates += 1
                errors.append((line_no, f"duplicate chromosome name {name!r} "
                                        f"(first seen on line {seen[name]})"))
            else:
                seen[name] = line_no

            try:
                size = int(size_token)
            except ValueError:
                errors.append((line_no, f"size is not an integer: {size_token!r}"))
                continue
            if size <= 0:
                errors.append((line_no, f"size must be >= 1, got {size}"))
                continue
            sizes.append(size)

    if records == 0:
        errors.append((0, "No chrom.sizes records found"))
    if order and order != sorted(order):
        warnings.append((0, "chromosomes are not in canonical order; the UCSC tools require the "
                            "coverage tracks to use the same order as this file"))
    if len(prefixes) > 1:
        warnings.append((0, "mixed chromosome naming ('chr'-prefixed and bare names); this breaks "
                            "track/annotation matching"))
    if space_separated:
        warnings.append((0, f"{space_separated} line(s) are not tab separated"))

    with open(report_path, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"input_compression\t{compression}\n")
        report.write(f"chromosomes\t{records}\n")
        report.write(f"total_bases\t{sum(sizes)}\n")
        report.write(f"min_size\t{min(sizes) if sizes else 0}\n")
        report.write(f"max_size\t{max(sizes) if sizes else 0}\n")
        report.write(f"duplicate_names\t{duplicates}\n")
        report.write(f"fai_shaped_lines\t{fai_columns}\n")
        report.write(f"overlong_names\t{long_names}\n")
        report.write(f"error_count\t{len(errors)}\n")
        report.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            report.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            report.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("chrom.sizes validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".sizes", ".chrom.sizes", ".txt")) as (in_fh, _), \
            open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
