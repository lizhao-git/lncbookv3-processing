"""Validate genePred / refFlat transcript tables (stdlib only)."""
from .io import open_text


def _detect_layout(cols):
    """Return 'genepred' or 'refflat' based on strand column position."""
    if len(cols) >= 3 and cols[2] in ("+", "-", "."):
        return "genepred"
    if len(cols) >= 4 and cols[3] in ("+", "-", "."):
        return "refflat"
    return None


def _parse_int_list(text, line_no, label, errors):
    values = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            values.append(int(token))
        except ValueError:
            errors.append((line_no, f"{label} contains a non-integer entry: {token}"))
    return values


def validate_genepred(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    records = 0
    chroms = []
    layouts = []
    compression = "plain"

    with open_text(input_path, preferred_exts=(".genepred", ".refflat")) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            records += 1
            cols = line.split("\t")
            layout = _detect_layout(cols)
            if layout is None:
                errors.append((line_no, "Cannot detect genePred/refFlat layout (strand column missing)"))
                continue
            if layout not in layouts:
                layouts.append(layout)
            if layout == "genepred":
                if len(cols) < 9:
                    errors.append((line_no, f"genePred requires >= 9 columns, got {len(cols)}"))
                    continue
                chrom, strand = cols[1], cols[2]
                tx_start, tx_end, cds_start, cds_end = cols[3], cols[4], cols[5], cols[6]
                exon_count, exon_starts, exon_ends = cols[7], cols[8], cols[9] if len(cols) > 9 else ""
            else:
                if len(cols) < 10:
                    errors.append((line_no, f"refFlat requires >= 10 columns, got {len(cols)}"))
                    continue
                chrom, strand = cols[2], cols[3]
                tx_start, tx_end, cds_start, cds_end = cols[4], cols[5], cols[6], cols[7]
                exon_count, exon_starts, exon_ends = cols[8], cols[9], cols[10] if len(cols) > 10 else ""

            if chrom not in chroms:
                chroms.append(chrom)
            if strand == ".":
                warnings.append((line_no, "Strand is unknown ('.')"))
            elif strand not in ("+", "-"):
                errors.append((line_no, f"Invalid strand: {strand}"))

            try:
                tx_start_i = int(tx_start)
                tx_end_i = int(tx_end)
                cds_start_i = int(cds_start)
                cds_end_i = int(cds_end)
            except ValueError:
                errors.append((line_no, "Transcript/CDS coordinates are not integers"))
                continue
            try:
                exon_count_i = int(exon_count)
            except ValueError:
                errors.append((line_no, "exonCount is not an integer"))
                continue

            if tx_start_i < 0:
                errors.append((line_no, "txStart must be >= 0"))
            if tx_end_i <= tx_start_i:
                errors.append((line_no, "txEnd must be > txStart"))
            if cds_start_i < tx_start_i or cds_end_i > tx_end_i:
                errors.append((line_no, "CDS bounds fall outside the transcript"))
            if cds_end_i < cds_start_i:
                errors.append((line_no, "cdsEnd must be >= cdsStart"))
            if exon_count_i < 1:
                errors.append((line_no, "exonCount must be >= 1"))

            starts = _parse_int_list(exon_starts, line_no, "exonStarts", errors)
            ends = _parse_int_list(exon_ends, line_no, "exonEnds", errors)
            if len(starts) != exon_count_i:
                errors.append((line_no, f"exonCount={exon_count_i} but found {len(starts)} exonStarts"))
            if len(ends) != exon_count_i:
                errors.append((line_no, f"exonCount={exon_count_i} but found {len(ends)} exonEnds"))
            if len(starts) == len(ends):
                for index, (exon_start, exon_end) in enumerate(zip(starts, ends)):
                    if exon_end <= exon_start:
                        errors.append((line_no, f"Exon {index + 1} end must be > start"))
                    if exon_start < tx_start_i or exon_end > tx_end_i:
                        errors.append((line_no, f"Exon {index + 1} falls outside the transcript bounds"))
                    if index and exon_start < starts[index - 1]:
                        warnings.append((line_no, f"Exon {index + 1} starts before exon {index} (unsorted)"))
            if exon_count_i >= 1 and starts and ends:
                if starts[0] != tx_start_i or ends[-1] != tx_end_i:
                    warnings.append((line_no, "Transcript bounds do not match the first/last exon"))

    if not layouts:
        warnings.append((0, "No records found"))
    elif len(layouts) > 1:
        warnings.append((0, f"Mixed layouts detected in one file: {', '.join(layouts)}"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"layout\t{layouts[0] if len(layouts) == 1 else 'mixed' if layouts else 'unknown'}\n")
        rep.write(f"records\t{records}\n")
        rep.write(f"chrom_count\t{len(chroms)}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("genePred/refFlat validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".genepred", ".refflat")) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
