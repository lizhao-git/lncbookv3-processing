"""Validate UCSC chain files used by liftOver (stdlib only)."""
from .io import open_text


def _parse_ints(tokens, line_no, label, errors):
    values = []
    for token in tokens:
        try:
            value = int(token)
        except ValueError:
            errors.append((line_no, f"{label} is not an integer: {token}"))
            return None
        if value < 0:
            errors.append((line_no, f"{label} must be >= 0: {token}"))
            return None
        values.append(value)
    return values


def validate_chain(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    chains = 0
    blocks = 0
    compression = "plain"
    current = None

    def close_chain(line_no, reason="blank line"):
        nonlocal current
        if current is None:
            return
        if current["open"]:
            errors.append(
                (line_no, f"Chain {current['id']} is not terminated by a size-only block line before {reason}")
            )
        if not current["open"]:
            t_span = current["t_used"]
            q_span = current["q_used"]
            if t_span != current["t_end"] - current["t_start"]:
                errors.append(
                    (line_no, f"Chain {current['id']} target blocks span {t_span} but header states "
                              f"{current['t_end'] - current['t_start']}")
                )
            if q_span != current["q_end"] - current["q_start"]:
                errors.append(
                    (line_no, f"Chain {current['id']} query blocks span {q_span} but header states "
                              f"{current['q_end'] - current['q_start']}")
                )
        current = None

    with open_text(input_path, preferred_exts=(".chain",)) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                close_chain(line_no)
                continue
            tokens = line.split()
            if tokens[0] == "chain":
                close_chain(line_no, reason="next chain header")
                if len(tokens) not in (12, 13):
                    errors.append((line_no, f"Chain header requires 13 fields (12 accepted legacy), got {len(tokens)}"))
                    current = None
                    continue
                if len(tokens) == 12:
                    warnings.append((line_no, "Chain header omits the trailing id field"))
                chains += 1
                score = _parse_ints([tokens[1]], line_no, "Score", errors)
                t_name, t_size, t_strand, t_start, t_end = tokens[2], tokens[3], tokens[4], tokens[5], tokens[6]
                q_name, q_size, q_strand, q_start, q_end = tokens[7], tokens[8], tokens[9], tokens[10], tokens[11]
                coords = _parse_ints([t_size, t_start, t_end, q_size, q_start, q_end], line_no, "Chain coordinate", errors)
                current = None
                if score is None or coords is None:
                    continue
                t_size_i, t_start_i, t_end_i, q_size_i, q_start_i, q_end_i = coords
                if t_strand not in ("+", "-"):
                    errors.append((line_no, f"Invalid target strand: {t_strand}"))
                if q_strand not in ("+", "-"):
                    errors.append((line_no, f"Invalid query strand: {q_strand}"))
                if t_size_i < 1 or q_size_i < 1:
                    errors.append((line_no, "Chain tSize/qSize must be >= 1"))
                if t_end_i <= t_start_i:
                    errors.append((line_no, "Chain tEnd must be > tStart"))
                if q_end_i <= q_start_i:
                    errors.append((line_no, "Chain qEnd must be > qStart"))
                if t_end_i > t_size_i:
                    errors.append((line_no, "Chain tEnd exceeds tSize"))
                if q_end_i > q_size_i:
                    errors.append((line_no, "Chain qEnd exceeds qSize"))
                current = {
                    "id": tokens[12] if len(tokens) == 13 else f"line{line_no}",
                    "t_start": t_start_i, "t_end": t_end_i, "t_used": 0,
                    "q_start": q_start_i, "q_end": q_end_i, "q_used": 0,
                    "open": True,
                    "names": f"{t_name}:{t_start_i}-{t_end_i} -> {q_name}:{q_start_i}-{q_end_i}",
                }
                continue

            if current is None:
                errors.append((line_no, "Block line appears before any chain header"))
                continue
            if not current["open"]:
                errors.append((line_no, f"Block line appears after chain {current['id']} was terminated"))
                continue
            values = _parse_ints(tokens, line_no, "Block size", errors)
            if values is None:
                continue
            if len(values) not in (1, 3):
                errors.append((line_no, f"Block line requires 1 or 3 integers, got {len(values)}"))
                continue
            blocks += 1
            if len(values) == 3:
                size, dt, dq = values
                current["t_used"] += size + dt
                current["q_used"] += size + dq
            else:
                current["t_used"] += values[0]
                current["q_used"] += values[0]
                current["open"] = False

    close_chain(0, reason="end of file")

    if chains == 0:
        errors.append((0, "No chain headers found"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"chains\t{chains}\n")
        rep.write(f"blocks\t{blocks}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("Chain validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".chain",)) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
