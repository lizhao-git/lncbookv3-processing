"""Minimal VCF reading, record model, writing, and validation helpers
(stdlib only).

Provides the shared plumbing for the VCF entry-point scripts: header
parsing, record iteration, INFO field handling, per-allele value
selection for multiallelic splits, genotype recoding, and file
validation.
"""
import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .io import open_text

_STRUCTURED_META_RE = re.compile(r"^##([A-Za-z][A-Za-z0-9_.]*)=<(.*)>\s*$")
_META_ATTR_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_.]*)=("[^"]*"|[^,]*)')

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def parse_structured_meta(line: str):
    """Parse ``##INFO=<...>``/``##FORMAT=<...>``/``##contig=<...>`` lines.

    Returns ``(kind, attributes)`` or None when the line is not a structured
    meta line. Quoted values are returned without their surrounding quotes.
    """
    match = _STRUCTURED_META_RE.match(line.strip())
    if not match:
        return None
    kind, body = match.groups()
    attrs = {}
    for key, value in _META_ATTR_RE.findall(body):
        attrs[key] = value.strip('"')
    return kind, attrs


def collect_field_definitions(meta_lines):
    """Return ``(info_defs, format_defs)``: field ID -> attribute dict."""
    info_defs = {}
    format_defs = {}
    for line in meta_lines:
        parsed = parse_structured_meta(line)
        if not parsed:
            continue
        kind, attrs = parsed
        name = attrs.get("ID")
        if not name:
            continue
        if kind == "INFO":
            info_defs[name] = attrs
        elif kind == "FORMAT":
            format_defs[name] = attrs
    return info_defs, format_defs


def allele_type(ref: str, alt: str) -> str:
    """Classify one ALT allele as snv, mnv, indel, or other."""
    if not ref or not alt or alt.startswith("<") or alt == "*" or "[" in alt or "]" in alt:
        return "other"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == len(alt):
        return "mnv"
    return "indel"


@dataclass
class VcfRecord:
    chrom: str
    pos: int
    id: str
    ref: str
    alts: List[str]
    qual: str
    filt: str
    info: str
    formats: List[str] = field(default_factory=list)
    samples: List[str] = field(default_factory=list)
    line_no: int = 0

    def alt_text(self) -> str:
        return ",".join(self.alts) if self.alts else "."

    def allele_types(self):
        """Set of allele types across all ALT alleles (``other`` when unset)."""
        if not self.alts:
            return {"other"}
        return {allele_type(self.ref, alt) for alt in self.alts}

    def to_line(self) -> str:
        cols = [self.chrom, str(self.pos), self.id, self.ref, self.alt_text(),
                self.qual, self.filt, self.info]
        if self.formats or self.samples:
            cols.append(":".join(self.formats))
            cols.extend(self.samples)
        return "\t".join(cols)


def parse_record(line: str, line_no: int = 0) -> VcfRecord:
    cols = line.split("\t")
    if len(cols) < 8:
        raise ValueError(f"expected at least 8 columns, got {len(cols)}")
    chrom, pos_text, id_, ref, alt_text, qual, filt, info = cols[:8]
    try:
        pos = int(pos_text)
    except ValueError:
        raise ValueError(f"POS is not an integer: {pos_text!r}")
    alts = [] if alt_text in (".", "") else alt_text.split(",")
    formats = []
    samples = []
    if len(cols) > 8:
        formats = [] if cols[8] in (".", "") else cols[8].split(":")
        samples = cols[9:]
    return VcfRecord(chrom, pos, id_, ref, alts, qual, filt, info, formats, samples, line_no)


def iter_vcf(path: str):
    """Yield ``(line_no, kind, item)`` tuples with kind 'meta'/'header'/'record'.

    Meta and header items are raw lines in file order; records are parsed
    :class:`VcfRecord` objects. Malformed records abort with a clear message.
    """
    with open_text(path) as (fh, _compression):
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                yield line_no, "meta", line
            elif line.startswith("#"):
                yield line_no, "header", line
            else:
                try:
                    record = parse_record(line, line_no)
                except ValueError as exc:
                    raise SystemExit(f"Invalid VCF record at line {line_no}: {exc}")
                yield line_no, "record", record


def header_samples(header_line: str) -> List[str]:
    cols = header_line.rstrip("\n").split("\t")
    if len(cols) <= 9:
        return []
    return cols[9:]


def build_header_line(samples) -> str:
    fixed = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if samples:
        fixed.append("FORMAT")
    return "\t".join(fixed + list(samples))


def parse_info(text: str):
    """Parse an INFO string into ``[(key, value_or_None)]``, order preserved."""
    if text in (".", ""):
        return []
    items = []
    for token in text.split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            items.append((key, value))
        else:
            items.append((token, None))
    return items


def format_info(items) -> str:
    if not items:
        return "."
    return ";".join(key if value is None else f"{key}={value}" for key, value in items)


def select_number_value(values: str, number: Optional[str], alt_index: int, alt_count: int):
    """Select the value(s) for one ALT allele of a multiallelic split.

    Follows the VCF Number conventions for ``A`` (one value per ALT) and
    ``R`` (REF + one per ALT); falls back to count-based inference when the
    header does not declare the field.
    """
    parts = values.split(",")
    if number == "A" and len(parts) >= alt_index + 1:
        return [parts[alt_index]]
    if number == "R" and len(parts) >= alt_index + 2:
        return [parts[0], parts[alt_index + 1]]
    if number in (None, "."):
        if len(parts) == alt_count:
            return [parts[alt_index]]
        if len(parts) == alt_count + 1:
            return [parts[0], parts[alt_index + 1]]
    return parts


def recode_gt(gt_text: str, alt_index: int) -> str:
    """Recode a GT value for a record split on ``alt_index`` (0-based).

    The kept ALT becomes allele ``1``; the REF stays ``0`` and every other
    ALT allele becomes missing (``.``). Separators are preserved.
    """
    out = []
    for token in re.split(r"([/|])", gt_text):
        if token in ("/", "|") or token == ".":
            out.append(token)
            continue
        if token.isdigit():
            value = int(token)
            if value == 0:
                out.append("0")
            elif value == alt_index + 1:
                out.append("1")
            else:
                out.append(".")
        else:
            out.append(".")
    return "".join(out)


def subset_genotype_field(values: str, alt_count: int, alt_index: int):
    """Return the diploid genotype-indexed subset for a biallelic split.

    ``values`` holds Number=G values ordered by the VCF convention
    (index = b*(b+1)/2 + a for a <= b). Returns None when the length does
    not match a diploid layout (caller should keep the field unchanged).
    """
    parts = values.split(",")
    expected = (alt_count + 1) * (alt_count + 2) // 2
    if len(parts) != expected:
        return None

    def index_of(a: int, b: int) -> int:
        if a > b:
            a, b = b, a
        return b * (b + 1) // 2 + a

    kept = alt_index + 1
    return [parts[index_of(0, 0)], parts[index_of(0, kept)], parts[index_of(kept, kept)]]


_VALID_BASES_RE = re.compile(r"^[ACGTNacgtn]+$")
_SYMBOLIC_ALT_RE = re.compile(r"^<([^<>]+)>$")
_BREAKEND_ALT_RE = re.compile(r"^[ACGTNacgtn.]*[\[\]][^\[\]]+[\[\]][ACGTNacgtn.]*$")
_VCF_HEADER_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


def validate_vcf(input_path: str, output_path: str, report_path: str):
    """Validate a VCF file and copy it through when no errors were found.

    Checks the ##fileformat / #CHROM headers, record structure (column
    counts, POS, REF, ALT incl. symbolic and breakend alleles, QUAL,
    FILTER, INFO) and per-sample FORMAT columns. Undeclared INFO/FORMAT
    IDs, FILTER ids, symbolic ALTs and contigs are reported as warnings,
    each warned name only once. Raises SystemExit when errors are found.
    """
    errors = []
    warnings = []
    has_fileformat = False
    has_header = False
    header_columns = 0
    samples = []
    info_ids = set()
    format_ids = set()
    alt_ids = set()
    filter_ids = set()
    contigs = set()
    warned_names = set()
    records = 0
    compression = "plain"

    def warn_once(line_no, name, message):
        key = (name, message)
        if key not in warned_names:
            warned_names.add(key)
            warnings.append((line_no, message))

    with open_text(input_path, preferred_exts=(".vcf",)) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                if line.startswith("##fileformat=VCF"):
                    has_fileformat = True
                parsed = parse_structured_meta(line)
                if parsed:
                    meta_kind, attrs = parsed
                    name = attrs.get("ID")
                    if not name:
                        continue
                    if meta_kind == "INFO":
                        info_ids.add(name)
                    elif meta_kind == "FORMAT":
                        format_ids.add(name)
                    elif meta_kind == "ALT":
                        alt_ids.add(name)
                    elif meta_kind == "FILTER":
                        filter_ids.add(name)
                    elif meta_kind == "contig":
                        contigs.add(name)
                continue
            if line.startswith("#"):
                if not line.startswith("#CHROM"):
                    warnings.append((line_no, "single-# comment line is not a #CHROM header"))
                    continue
                if has_header:
                    errors.append((line_no, "Duplicate #CHROM header line"))
                has_header = True
                cols = line.split("\t")
                if len(cols) < 8:
                    errors.append((line_no, f"#CHROM header has fewer than 8 columns: {len(cols)}"))
                    continue
                if cols[:8] != _VCF_HEADER_COLUMNS:
                    errors.append((line_no, "#CHROM header columns do not match the VCF specification"))
                header_columns = len(cols)
                if len(cols) > 8:
                    if cols[8] != "FORMAT":
                        errors.append((line_no, "Sample columns are present but column 9 is not FORMAT"))
                    samples = cols[9:]
                continue

            if not has_header:
                errors.append((line_no, "Variant record before the #CHROM header line"))
            records += 1
            cols = line.split("\t")
            if len(cols) < 8:
                errors.append((line_no, f"Expected at least 8 columns, got {len(cols)}"))
                continue
            if header_columns and len(cols) != header_columns:
                errors.append(
                    (line_no, f"Record has {len(cols)} columns but the header declares {header_columns}")
                )
            chrom, pos_text, _id, ref, alt_text, qual, filt, info = cols[:8]

            if contigs and chrom not in contigs:
                warn_once(line_no, chrom, f"CHROM not declared in ##contig header: {chrom}")
            try:
                if int(pos_text) <= 0:
                    errors.append((line_no, "POS must be > 0"))
            except ValueError:
                errors.append((line_no, f"POS is not an integer: {pos_text!r}"))

            if not _VALID_BASES_RE.match(ref):
                errors.append((line_no, f"REF must be a non-empty sequence of A/C/G/T/N bases: {ref!r}"))
            for token in ([] if alt_text in (".", "") else alt_text.split(",")):
                if token in ("*",) or _VALID_BASES_RE.match(token):
                    continue
                symbolic = _SYMBOLIC_ALT_RE.match(token)
                if symbolic:
                    alt_id = symbolic.group(1)
                    warn_once(line_no, alt_id,
                              f"Symbolic ALT not declared in the ##ALT header: <{alt_id}>")
                    continue
                if _BREAKEND_ALT_RE.match(token):
                    continue
                errors.append((line_no, f"Invalid ALT allele: {token!r}"))
            if alt_text in (".", ""):
                warnings.append((line_no, "ALT is '.' (no variant allele)"))

            if qual != ".":
                try:
                    if float(qual) < 0:
                        errors.append((line_no, "QUAL must be >= 0"))
                except ValueError:
                    errors.append((line_no, f"QUAL is not numeric: {qual!r}"))

            if filt in ("", "."):
                pass
            else:
                for filt_id in filt.split(";"):
                    if filt_id != "PASS" and filt_id not in filter_ids:
                        warn_once(line_no, filt_id, f"FILTER id not declared in header: {filt_id}")

            for key, _value in parse_info(info):
                if not key:
                    errors.append((line_no, "INFO field with an empty key"))
                elif key not in info_ids:
                    warn_once(line_no, key, f"INFO key not declared in header: {key}")

            if len(cols) > 8:
                format_text = cols[8]
                fmt_fields = [] if format_text in (".", "") else format_text.split(":")
                if fmt_fields:
                    if "GT" in fmt_fields and fmt_fields[0] != "GT":
                        warnings.append((line_no, "GT is present but is not the first FORMAT field"))
                    for field in fmt_fields:
                        if field not in format_ids:
                            warn_once(line_no, field, f"FORMAT field not declared in header: {field}")

    if not has_fileformat:
        errors.append((0, "Missing ##fileformat=VCF header"))
    if not has_header:
        errors.append((0, "Missing #CHROM header line"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"has_fileformat_header\t{has_fileformat}\n")
        rep.write(f"has_chrom_header\t{has_header}\n")
        rep.write(f"samples\t{len(samples)}\n")
        rep.write(f"contigs\t{len(contigs)}\n")
        rep.write(f"info_fields\t{len(info_ids)}\n")
        rep.write(f"format_fields\t{len(format_ids)}\n")
        rep.write(f"records\t{records}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        rep.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            rep.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("VCF validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".vcf",)) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
