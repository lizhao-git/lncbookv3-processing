"""Validate PSL alignment files produced by BLAT and the UCSC tools (stdlib only).

The structural scan below runs everywhere; when the kent ``pslCheck`` binary is
available it is run as well and its verdict is folded into the same report.
``pslCheck`` aborts on a record whose block list length disagrees with
``blockCount`` (``pslLoad`` assertion), so the structural pass runs first and
suppresses the call when a record cannot be parsed safely.
"""
import os
import re
import shutil
import subprocess

from .io import open_text

PSL_FIELD_COUNT = 21
PSLX_FIELD_COUNT = 23
STRANDS = ("+", "-", "++", "+-", "--", "-+")

HEADER_KEYWORDS = ("psLayout", "match", "score")

INTEGER_FIELDS = (
    (0, "matches"),
    (1, "misMatches"),
    (2, "repMatches"),
    (3, "nCount"),
    (4, "qNumInsert"),
    (5, "qBaseInsert"),
    (6, "tNumInsert"),
    (7, "tBaseInsert"),
    (10, "qSize"),
    (11, "qStart"),
    (12, "qEnd"),
    (14, "tSize"),
    (15, "tStart"),
    (16, "tEnd"),
    (17, "blockCount"),
)
LIST_FIELDS = ((18, "blockSizes"), (19, "qStarts"), (20, "tStarts"))
NAME_FIELDS = ((9, "qName"), (13, "tName"))

PSLCHECK_SUMMARY = re.compile(r"checked:\s*(\d+)\s+failed:\s*(\d+)\s+errors:\s*(\d+)")


def _parse_int(token, line_no, label, errors):
    try:
        value = int(token)
    except ValueError:
        errors.append((line_no, f"{label} is not an integer: {token!r}"))
        return None
    if value < 0:
        errors.append((line_no, f"{label} must be >= 0: {token!r}"))
        return None
    return value


def _parse_int_list(token, line_no, label, errors):
    """Parse a BLAT comma-terminated list such as ``10,20,``."""
    values = []
    for part in token.split(","):
        if part == "":
            continue
        try:
            value = int(part)
        except ValueError:
            errors.append((line_no, f"{label} has a non-integer element: {part!r}"))
            return None
        if value < 0:
            errors.append((line_no, f"{label} has a negative element: {part!r}"))
            return None
        values.append(value)
    if not values:
        errors.append((line_no, f"{label} is empty: {token!r}"))
        return None
    return values


def _is_header_line(line, tokens):
    """True for the psLayout label rows and the dashed rule."""
    if line.replace("-", "").strip() == "":
        return True
    return bool(tokens) and tokens[0] in HEADER_KEYWORDS


def _check_record(fields, line_no, errors, warnings, state):
    """Field-level checks for one PSL record. Returns a key for duplicate detection."""
    numeric = {}
    parseable = True
    for index, label in INTEGER_FIELDS:
        value = _parse_int(fields[index], line_no, label, errors)
        if value is None:
            parseable = False
        else:
            numeric[label] = value

    lists = {}
    for index, label in LIST_FIELDS:
        if not fields[index].endswith(","):
            state["missing_trailing_comma"] += 1
        values = _parse_int_list(fields[index], line_no, label, errors)
        if values is None:
            parseable = False
        else:
            lists[label] = values

    strand = fields[8]
    if strand not in STRANDS:
        errors.append((line_no, f"strand {strand!r} is not one of {', '.join(STRANDS)}"))
    state["strands"][strand] = state["strands"].get(strand, 0) + 1

    for index, label in NAME_FIELDS:
        if not fields[index]:
            errors.append((line_no, f"{label} is empty"))

    if not parseable:
        state["pslcheck_unsafe"] = True
        return None

    block_count = numeric["blockCount"]
    block_sizes = lists["blockSizes"]
    q_starts = lists["qStarts"]
    t_starts = lists["tStarts"]
    state["blocks"] += sum(block_sizes)

    if block_count < 1:
        errors.append((line_no, "blockCount must be >= 1"))
        state["pslcheck_unsafe"] = True
        return None
    if not (len(block_sizes) == len(q_starts) == len(t_starts) == block_count):
        errors.append(
            (line_no, f"blockCount is {block_count} but blockSizes/qStarts/tStarts hold "
                      f"{len(block_sizes)}/{len(q_starts)}/{len(t_starts)} entries")
        )
        # pslLoad asserts on this and aborts, so pslCheck must not be called.
        state["pslcheck_unsafe"] = True
        return None

    q_size, q_start, q_end = numeric["qSize"], numeric["qStart"], numeric["qEnd"]
    t_size, t_start, t_end = numeric["tSize"], numeric["tStart"], numeric["tEnd"]

    if q_size < 1:
        errors.append((line_no, "qSize must be >= 1"))
    if t_size < 1:
        errors.append((line_no, "tSize must be >= 1"))
    if q_end < q_start:
        errors.append((line_no, f"qEnd ({q_end}) must be >= qStart ({q_start})"))
    if t_end < t_start:
        errors.append((line_no, f"tEnd ({t_end}) must be >= tStart ({t_start})"))
    if q_end > q_size:
        errors.append((line_no, f"qEnd ({q_end}) exceeds qSize ({q_size})"))
    if t_end > t_size:
        errors.append((line_no, f"tEnd ({t_end}) exceeds tSize ({t_size})"))

    if any(q_starts[i] > q_starts[i + 1] for i in range(block_count - 1)):
        errors.append((line_no, "qStarts are not in ascending order"))
    if any(t_starts[i] > t_starts[i + 1] for i in range(block_count - 1)):
        errors.append((line_no, "tStarts are not in ascending order"))

    for i, (size, qs, ts) in enumerate(zip(block_sizes, q_starts, t_starts)):
        if qs < q_start:
            errors.append((line_no, f"block {i} qStart ({qs}) precedes the record qStart ({q_start})"))
        if ts < t_start:
            errors.append((line_no, f"block {i} tStart ({ts}) precedes the record tStart ({t_start})"))
        if qs + size > q_size:
            errors.append((line_no, f"block {i} qStart+blockSize ({qs + size}) exceeds qSize ({q_size})"))
        if ts + size > t_size:
            errors.append((line_no, f"block {i} tStart+blockSize ({ts + size}) exceeds tSize ({t_size})"))

    if q_starts[0] != q_start:
        warnings.append((line_no, f"qStart ({q_start}) differs from the first qStarts entry ({q_starts[0]})"))
    if q_starts[-1] + block_sizes[-1] != q_end:
        warnings.append((line_no, f"qEnd ({q_end}) differs from the last block end "
                                  f"({q_starts[-1] + block_sizes[-1]})"))

    aligned = sum(block_sizes)
    counted = numeric["matches"] + numeric["misMatches"] + numeric["repMatches"] + numeric["nCount"]
    if counted != aligned:
        message = (f"matches+misMatches+repMatches+nCount ({counted}) does not match the alignment "
                   f"size ({aligned})")
        # The identity is exact for nucleotide PSLs; translated PSLs report protein counts.
        if len(strand) == 1:
            errors.append((line_no, message))
        else:
            warnings.append((line_no, message + "; expected for translated PSL records"))

    # The unaligned bases and the gap counts are both derivable from the block
    # layout, so all four Q/T gap fields can be checked exactly rather than
    # compared against a blockCount-1 rule of thumb.
    for side, span_start, span_end, starts in (
        ("q", q_start, q_end, q_starts),
        ("t", t_start, t_end, t_starts),
    ):
        leading = starts[0] - span_start
        trailing = span_end - (starts[-1] + block_sizes[-1])
        internal = [starts[i + 1] - (starts[i] + block_sizes[i]) for i in range(block_count - 1)]
        gap_count = (1 if leading else 0) + sum(1 for gap in internal if gap) + (1 if trailing else 0)

        base_insert = numeric[f"{side}BaseInsert"]
        if base_insert != leading + trailing + sum(internal):
            errors.append((line_no, f"{side}BaseInsert ({base_insert}) does not match the unaligned "
                                    f"bases ({leading + trailing + sum(internal)})"))
        num_insert = numeric[f"{side}NumInsert"]
        if num_insert != gap_count:
            message = (f"{side}NumInsert ({num_insert}) does not match the {gap_count} {side} gap(s) "
                       f"implied by the block layout")
            # BLAT computes these counts inconsistently for translated PSL records.
            if len(strand) == 1:
                errors.append((line_no, message))
            else:
                warnings.append((line_no, message + "; expected for translated PSL records"))

    state["queries"].add(fields[9])
    state["targets"].add(fields[13])
    return (fields[9], fields[13], q_start, t_start, block_count, tuple(block_sizes))


def _run_pslcheck(input_path, report_dir, target_sizes, query_sizes):
    """Run the kent pslCheck splitter. Returns None when the binary is absent."""
    tool = shutil.which("pslCheck")
    if tool is None:
        return None

    pass_path = os.path.join(report_dir, "pslcheck_pass.psl")
    fail_path = os.path.join(report_dir, "pslcheck_fail.psl")
    command = [tool, "-filter", f"-pass={pass_path}", f"-fail={fail_path}"]
    if target_sizes:
        command.append(f"-targetSizes={target_sizes}")
    if query_sizes:
        command.append(f"-querySizes={query_sizes}")
    command.append(input_path)

    result = subprocess.run(command, capture_output=True, text=True)
    # pslCheck reports both its per-record errors and the trailing "checked:"
    # summary on stderr, so the two have to be told apart.
    stream = f"{result.stdout or ''}\n{result.stderr or ''}"
    summary = PSLCHECK_SUMMARY.search(stream)
    detail = []
    for line in (result.stderr or "").splitlines():
        stripped = line.strip()
        if not stripped or PSLCHECK_SUMMARY.search(stripped):
            continue
        if stripped.startswith("Assertion failed"):
            continue
        detail.append(stripped)
    outcome = {
        "returncode": result.returncode,
        "checked": int(summary.group(1)) if summary else None,
        "failed": int(summary.group(2)) if summary else None,
        "errors": int(summary.group(3)) if summary else None,
        "detail": detail,
    }
    for path in (pass_path, fail_path):
        try:
            os.remove(path)
        except OSError:
            pass
    return outcome


def validate_psl(input_path, output_path, report_path, target_sizes=None, query_sizes=None,
                 use_pslcheck=True):
    errors = []
    warnings = []
    compression = "plain"
    state = {
        "pslcheck_unsafe": False,
        "blocks": 0,
        "missing_trailing_comma": 0,
        "queries": set(),
        "targets": set(),
        "strands": {},
    }
    header_lines = 0
    records = 0
    pslx_records = 0
    space_separated = 0
    duplicates = 0
    seen = set()

    with open_text(input_path, preferred_exts=(".psl",)) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            if _is_header_line(line, tokens):
                if records:
                    # .cat-style repeats are tolerated; a stray rule after the data is not.
                    warnings.append((line_no, "psLayout header line found after PSL records"))
                header_lines += 1
                continue
            if not tokens[0].isdigit():
                errors.append((line_no, f"PSL record must start with an integer match count: {tokens[0]!r}"))
                state["pslcheck_unsafe"] = True
                continue

            fields = line.split("\t")
            if len(fields) != len(tokens):
                space_separated += 1
                fields = tokens
            if len(fields) not in (PSL_FIELD_COUNT, PSLX_FIELD_COUNT):
                errors.append((line_no, f"PSL record needs {PSL_FIELD_COUNT} fields "
                                        f"({PSLX_FIELD_COUNT} for psLx), got {len(fields)}"))
                state["pslcheck_unsafe"] = True
                continue
            if len(fields) == PSL_FIELD_COUNT + 1:
                warnings.append((line_no, "PSL record has 22 fields; expected 21 (PSL) or 23 (psLx)"))

            records += 1
            if len(fields) == PSLX_FIELD_COUNT:
                pslx_records += 1
                for index, label in ((21, "qStrand"), (22, "tStrand")):
                    if fields[index] not in ("+", "-"):
                        errors.append((line_no, f"{label} must be '+' or '-', got {fields[index]!r}"))

            key = _check_record(fields, line_no, errors, warnings, state)
            if key is not None:
                if key in seen:
                    duplicates += 1
                seen.add(key)

    if records == 0:
        errors.append((0, "No PSL records found"))
    if not header_lines:
        warnings.append((0, "No psLayout header block; treating the file as headerless PSL"))
    if space_separated:
        warnings.append((0, f"{space_separated} record(s) are not tab separated"))

    pslcheck = None
    if use_pslcheck and not state["pslcheck_unsafe"]:
        pslcheck = _run_pslcheck(
            input_path, os.path.dirname(os.path.abspath(report_path)),
            target_sizes, query_sizes,
        )
        if pslcheck is not None:
            for line in pslcheck["detail"]:
                errors.append((0, line))
            if pslcheck["returncode"] != 0:
                errors.append((0, f"pslCheck exited with status {pslcheck['returncode']}"))
            elif pslcheck["failed"]:
                errors.append((0, f"pslCheck rejected {pslcheck['failed']} record(s)"))
    skipped_reason = ""
    if not use_pslcheck:
        skipped_reason = "disabled with --skip-pslcheck"
    elif state["pslcheck_unsafe"]:
        skipped_reason = "a record could not be parsed safely (pslCheck would abort)"

    with open(report_path, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"input_compression\t{compression}\n")
        report.write(f"header_lines\t{header_lines}\n")
        report.write(f"records\t{records}\n")
        report.write(f"pslx_records\t{pslx_records}\n")
        report.write(f"query_count\t{len(state['queries'])}\n")
        report.write(f"target_count\t{len(state['targets'])}\n")
        report.write(f"alignment_blocks\t{state['blocks']}\n")
        report.write("strand_counts\t" + " | ".join(
            f"{s}:{n}" for s, n in sorted(state["strands"].items())) + "\n")
        report.write(f"duplicate_records\t{duplicates}\n")
        report.write(f"records_missing_trailing_comma\t{state['missing_trailing_comma']}\n")
        if pslcheck is not None:
            report.write("pslcheck\tran\n")
            report.write(f"pslcheck_checked\t{pslcheck['checked']}\n")
            report.write(f"pslcheck_failed\t{pslcheck['failed']}\n")
        else:
            report.write("pslcheck\tnot_run\n")
            report.write(f"pslcheck_skipped_reason\t{skipped_reason}\n")
        report.write(f"error_count\t{len(errors)}\n")
        report.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            report.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            report.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("PSL validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".psl",)) as (in_fh, _), \
            open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
