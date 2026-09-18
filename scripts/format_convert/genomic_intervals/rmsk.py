"""Validate RepeatMasker .out annotation files (stdlib only).

Handles the plain ``.out`` layout, the ``-a`` alignment variant (each hit is
followed by a continuation line repeating the coordinates) and the ``.cat``
layout, where the two title lines are repeated before every sequence.
"""
import re

from .io import open_text

TITLE_KEYWORDS = ("SW", "score")
TERMINATOR = "*"
NO_HITS_MESSAGE = "There were no repetitive sequences detected"

PRIMARY_FIELDS = (14, 15)
CONTINUATION_FIELDS = (10, 11)
STRANDS = ("+", "-", "C", "c")

LEFT_FIELD = re.compile(r"^\((-?\d+)\)$")
HEADER_FIELD = 0


def _parse_int(token, line_no, label, errors, minimum=None):
    try:
        value = int(token)
    except ValueError:
        errors.append((line_no, f"{label} is not an integer: {token!r}"))
        return None
    if minimum is not None and value < minimum:
        errors.append((line_no, f"{label} must be >= {minimum}, got {value}"))
        return None
    return value


def _parse_percent(token, line_no, label, errors):
    try:
        value = float(token)
    except ValueError:
        errors.append((line_no, f"{label} is not numeric: {token!r}"))
        return None
    if not 0.0 <= value <= 100.0:
        errors.append((line_no, f"{label} must be between 0 and 100, got {token!r}"))
        return None
    return value


def _parse_left(token, line_no, label, errors, allow_negative=False):
    """Parse a RepeatMasker '(N)' remaining-bases field."""
    match = LEFT_FIELD.match(token)
    if not match:
        errors.append((line_no, f"{label} must look like '(123)', got {token!r}"))
        return None
    value = int(match.group(1))
    if value < 0 and not allow_negative:
        errors.append((line_no, f"{label} must be >= 0, got {token!r}"))
        return None
    return value


def _check_halves(fields, offset, line_no, errors, state):
    """Validate the repeat half (name, class, coordinates, ID) of a record."""
    repeat, repeat_class = fields[offset], fields[offset + 1]
    if not repeat:
        errors.append((line_no, "repeat name is empty"))
    if not repeat_class:
        errors.append((line_no, "repeat class/family is empty"))
    r_begin = _parse_int(fields[offset + 2], line_no, "repeat begin", errors, minimum=1)
    r_end = _parse_int(fields[offset + 3], line_no, "repeat end", errors, minimum=1)
    r_left = _parse_left(fields[offset + 4], line_no, "repeat (left)", errors, allow_negative=True)
    record_id = _parse_int(fields[offset + 5], line_no, "ID", errors, minimum=1)

    if r_begin is not None and r_end is not None and r_end < r_begin:
        errors.append((line_no, f"repeat end ({r_end}) must be >= repeat begin ({r_begin})"))

    state["repeats"].add(repeat)
    if repeat_class:
        state["classes"].add(repeat_class)
    if r_end is not None and r_left is not None:
        # (left) is measured against the consensus, so the sum is the consensus
        # length; repeats that overrun the consensus report a negative (left).
        consensus = r_end + r_left
        known = state["repeat_lengths"].get(repeat)
        if known is None:
            state["repeat_lengths"][repeat] = consensus
        elif known != consensus and r_left >= 0:
            state["repeat_length_mismatch"] += 1
    if record_id is not None:
        if state["last_id"] is None:
            state["first_id"] = record_id
        elif record_id <= state["last_id"]:
            state["id_order_breaks"] += 1
        state["last_id"] = record_id
    return record_id


def _register_query(name, begin, end, left, line_no, errors, state):
    """Track query coordinates and verify qEnd + (left) equals the query length."""
    if not name:
        errors.append((line_no, "query sequence name is empty"))
    if end is not None and left is not None:
        length = end + left
        known = state["query_lengths"].get(name)
        if known is None:
            state["query_lengths"][name] = length
        elif known != length:
            errors.append((line_no, f"query {name!r} implies length {length} but earlier records "
                                    f"imply {known} (end + (left) must be constant)"))
    if begin is not None and end is not None:
        previous_end = state["query_last_end"].get(name)
        if previous_end is not None:
            if begin < state["query_last_begin"].get(name, begin):
                state["query_order_breaks"] += 1
            if begin <= previous_end:
                state["query_overlaps"] += 1
        state["query_last_end"][name] = end
        state["query_last_begin"][name] = begin
    state["queries"].add(name)
    return name


def _check_primary(fields, line_no, errors, warnings, state):
    """Validate one summary record (leading SW score + four percentage columns)."""
    has_strand = len(fields) == 15
    expected = 15 if has_strand else 14
    if len(fields) not in PRIMARY_FIELDS:
        errors.append((line_no, f"RepeatMasker record needs {expected} fields, got {len(fields)}"))
        return None

    _parse_int(fields[0], line_no, "SW score", errors, minimum=0)
    for index, label in ((1, "perc div."), (2, "perc del."), (3, "perc ins.")):
        _parse_percent(fields[index], line_no, label, errors)

    name = fields[4]
    begin = _parse_int(fields[5], line_no, "query begin", errors, minimum=1)
    end = _parse_int(fields[6], line_no, "query end", errors, minimum=1)
    left = _parse_left(fields[7], line_no, "query (left)", errors)
    if begin is not None and end is not None and end < begin:
        errors.append((line_no, f"query end ({end}) must be >= query begin ({begin})"))

    if has_strand:
        strand = fields[8]
        if strand not in STRANDS:
            errors.append((line_no, f"strand {strand!r} is not one of {', '.join(STRANDS)}"))
        state["strands"][strand] = state["strands"].get(strand, 0) + 1
        offset = 9
    else:
        state["strands"]["+"] = state["strands"].get("+", 0) + 1
        offset = 8

    _register_query(name, begin, end, left, line_no, errors, state)
    record_id = _check_halves(fields, offset, line_no, errors, state)

    state["last_primary"] = {
        "name": name, "begin": begin, "end": end, "left": left,
        "strand": fields[8] if has_strand else "+",
        "repeat": fields[offset], "repeat_class": fields[offset + 1],
        "r_begin": fields[offset + 2], "r_end": fields[offset + 3],
        "r_left": fields[offset + 4], "id": fields[offset + 5],
        "line": line_no,
    }
    return record_id


def _check_continuation(fields, line_no, errors, state):
    """Validate an '-a' alignment line, which must restate the preceding record."""
    last = state["last_primary"]
    if last is None:
        errors.append((line_no, "alignment continuation line appears before any summary record"))
        return
    if len(fields) not in CONTINUATION_FIELDS:
        errors.append((line_no, f"alignment continuation line needs 10 fields, got {len(fields)}"))
        return

    has_strand = len(fields) == 11
    offset = 5 if has_strand else 4
    expected = (
        ("query name", fields[0], last["name"]),
        ("query begin", fields[1], str(last["begin"])),
        ("query end", fields[2], str(last["end"])),
        ("query (left)", fields[3], f"({last['left']})"),
        ("repeat name", fields[offset], last["repeat"]),
        ("repeat begin", fields[offset + 2], str(last["r_begin"])),
        ("repeat end", fields[offset + 3], str(last["r_end"])),
        ("repeat (left)", fields[offset + 4], last["r_left"]),
        ("ID", fields[offset + 5], str(last["id"])),
    )
    for label, actual, wanted in expected:
        if actual != wanted:
            errors.append((line_no, f"continuation {label} is {actual!r} but line {last['line']} "
                                    f"states {wanted!r}"))


def validate_rmsk(input_path, output_path, report_path):
    errors = []
    warnings = []
    compression = "plain"
    state = {
        "queries": set(),
        "repeats": set(),
        "classes": set(),
        "query_lengths": {},
        "query_last_end": {},
        "query_last_begin": {},
        "repeat_lengths": {},
        "repeat_length_mismatch": 0,
        "query_order_breaks": 0,
        "query_overlaps": 0,
        "id_order_breaks": 0,
        "last_id": None,
        "first_id": "",
        "last_primary": None,
        "strands": {},
    }
    header_lines = 0
    preamble_lines = 0
    records = 0
    continuation_lines = 0
    terminators = 0
    no_hits_messages = 0

    with open_text(input_path, preferred_exts=(".out", ".cat")) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line == TERMINATOR:
                terminators += 1
                continue
            tokens = line.split()
            if tokens[HEADER_FIELD] in TITLE_KEYWORDS:
                header_lines += 1
                continue
            if line.startswith(NO_HITS_MESSAGE):
                no_hits_messages += 1
                continue

            if tokens[HEADER_FIELD].isdigit():
                records += 1
                _check_primary(tokens, line_no, errors, warnings, state)
                continue

            if records == 0 and len(tokens) not in CONTINUATION_FIELDS:
                preamble_lines += 1
                continue
            continuation_lines += 1
            _check_continuation(tokens, line_no, errors, state)

    if records == 0:
        if no_hits_messages:
            warnings.append((0, "No repeat records; the file reports no repetitive sequences"))
        else:
            errors.append((0, "No RepeatMasker records found"))
    if preamble_lines:
        warnings.append((0, f"{preamble_lines} unrecognised line(s) before the first record were skipped"))
    if state["query_order_breaks"]:
        warnings.append((0, f"{state['query_order_breaks']} record(s) are out of order within their query"))
    if state["query_overlaps"]:
        warnings.append((0, f"{state['query_overlaps']} record(s) overlap the previous hit on the same query"))
    if state["id_order_breaks"]:
        warnings.append((0, f"{state['id_order_breaks']} record ID(s) are not strictly increasing"))
    if state["repeat_length_mismatch"]:
        warnings.append((0, f"{state['repeat_length_mismatch']} record(s) disagree with the consensus "
                            f"length implied for the same repeat"))

    with open(report_path, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"input_compression\t{compression}\n")
        report.write(f"header_lines\t{header_lines}\n")
        report.write(f"records\t{records}\n")
        report.write(f"continuation_lines\t{continuation_lines}\n")
        report.write(f"terminators\t{terminators}\n")
        report.write(f"no_hits_messages\t{no_hits_messages}\n")
        report.write(f"query_count\t{len(state['queries'])}\n")
        report.write(f"repeat_count\t{len(state['repeats'])}\n")
        report.write(f"class_count\t{len(state['classes'])}\n")
        report.write(f"first_id\t{state['first_id']}\n")
        report.write(f"last_id\t{state['last_id'] if state['last_id'] is not None else ''}\n")
        report.write("strand_counts\t" + " | ".join(
            f"{s}:{n}" for s, n in sorted(state["strands"].items())) + "\n")
        report.write(f"error_count\t{len(errors)}\n")
        report.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            report.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            report.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("RepeatMasker validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".out", ".cat")) as (in_fh, _), \
            open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
