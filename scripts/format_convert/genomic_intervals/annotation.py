import re
from collections import Counter, defaultdict
from dataclasses import dataclass

from .io import open_text


GTF_EXTS = (".gtf", ".gff", ".gff3")


@dataclass
class AnnotationRecord:
    chrom: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attrs_text: str
    attrs: dict


def detect_annotation_format(path: str, explicit_format: str = "auto") -> str:
    if explicit_format != "auto":
        return explicit_format
    lower = path.lower()
    if lower.endswith((".gff", ".gff3")):
        return "gff3"
    return "gtf"


def parse_gtf_attrs(attr_text: str) -> dict:
    attrs = {}
    for match in re.finditer(r'([^\s;]+)\s+"([^"]*)"\s*;?', attr_text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def parse_gff3_attrs(attr_text: str) -> dict:
    attrs = {}
    for token in attr_text.split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            attrs[key] = value
        else:
            attrs[token] = ""
    return attrs


def parse_attrs(attr_text: str, fmt: str) -> dict:
    if fmt == "gff3":
        return parse_gff3_attrs(attr_text)
    return parse_gtf_attrs(attr_text)


_GFF3_ESCAPE_RE = re.compile(r"[%;=,&\t\n\r]")


def gff3_escape(value: str) -> str:
    """Percent-encode GFF3-reserved characters in an attribute value."""
    return _GFF3_ESCAPE_RE.sub(lambda match: f"%{ord(match.group(0)):02X}", value)


def gff3_unescape(value: str) -> str:
    """Decode percent-encoded sequences in a GFF3 attribute value."""
    return re.sub(r"%([0-9A-Fa-f]{2})", lambda match: chr(int(match.group(1), 16)), value)


def iter_annotation(path: str, fmt: str = "auto"):
    fmt = detect_annotation_format(path, fmt)
    with open_text(path, preferred_exts=GTF_EXTS) as (fh, _compression):
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                yield line_no, None, f"Expected 9 columns, got {len(cols)}"
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs_text = cols
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                yield line_no, None, "Start/End are not integers"
                continue
            yield line_no, AnnotationRecord(
                chrom=chrom,
                source=source,
                feature=feature,
                start=start_i,
                end=end_i,
                score=score,
                strand=strand,
                frame=frame,
                attrs_text=attrs_text,
                attrs=parse_attrs(attrs_text, fmt),
            ), None


def feature_id(record: AnnotationRecord) -> str:
    attrs = record.attrs
    return (
        attrs.get("gene_id")
        or attrs.get("transcript_id")
        or attrs.get("ID")
        or attrs.get("Name")
        or "NA"
    )


def gene_id(record: AnnotationRecord) -> str:
    attrs = record.attrs
    return attrs.get("gene_id") or attrs.get("gene") or attrs.get("Parent") or attrs.get("ID") or "NA"


def transcript_id(record: AnnotationRecord) -> str:
    attrs = record.attrs
    return attrs.get("transcript_id") or attrs.get("transcript") or attrs.get("ID") or attrs.get("Parent") or "NA"


def gene_type(record: AnnotationRecord) -> str:
    attrs = record.attrs
    return attrs.get("gene_type") or attrs.get("gene_biotype") or attrs.get("biotype") or "NA"


def transcript_type(record: AnnotationRecord) -> str:
    attrs = record.attrs
    return attrs.get("transcript_type") or attrs.get("transcript_biotype") or attrs.get("biotype") or "NA"


def to_bed_row(record: AnnotationRecord, kind: str, fid: str, gid: str, tid: str, gtype: str, ttype: str):
    return [
        record.chrom,
        str(max(0, record.start - 1)),
        str(record.end),
        kind,
        fid,
        gid,
        tid,
        record.strand,
        gtype,
        ttype,
    ]


def validate_annotation(input_path: str, output_path: str, report_path: str, fmt: str = "auto"):
    fmt = detect_annotation_format(input_path, fmt)
    errors = []
    feature_counter = Counter()
    line_count = 0
    compression = "plain"

    with open_text(input_path, preferred_exts=GTF_EXTS) as (fh, detected_compression):
        compression = detected_compression
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            line_count += 1
            cols = line.split("\t")
            if len(cols) != 9:
                errors.append((line_no, f"Expected 9 columns, got {len(cols)}"))
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs_text = cols
            feature_counter[feature] += 1
            try:
                start_i = int(start)
                end_i = int(end)
                if start_i <= 0 or end_i <= 0:
                    errors.append((line_no, "Start/End must be positive integers"))
                if start_i > end_i:
                    errors.append((line_no, "Start must be <= End"))
            except ValueError:
                errors.append((line_no, "Start/End are not integers"))
            if strand not in {"+", "-", ".", "?"}:
                errors.append((line_no, f"Invalid strand value: {strand}"))
            attrs = parse_attrs(attrs_text, fmt)
            if feature == "gene" and not (attrs.get("gene_id") or attrs.get("ID")):
                errors.append((line_no, "gene feature missing gene_id/ID"))
            if feature in {"transcript", "mRNA", "lnc_RNA", "exon", "intron"}:
                if not (attrs.get("gene_id") or attrs.get("Parent") or attrs.get("gene")):
                    errors.append((line_no, f"{feature} missing gene_id/Parent"))
                if feature in {"transcript", "mRNA", "lnc_RNA", "exon", "intron"} and not (
                    attrs.get("transcript_id") or attrs.get("ID") or attrs.get("Parent")
                ):
                    errors.append((line_no, f"{feature} missing transcript_id/ID/Parent"))

    with open(report_path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"format\t{fmt}\n")
        rep.write(f"input_compression\t{compression}\n")
        rep.write(f"parsed_records\t{line_count}\n")
        for name, count in sorted(feature_counter.items()):
            rep.write(f"feature_{name}\t{count}\n")
        rep.write(f"error_count\t{len(errors)}\n")
        if errors:
            rep.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")

    if errors:
        raise SystemExit("Annotation validation failed. See report for details.")

    with open_text(input_path, preferred_exts=GTF_EXTS) as (in_fh, _), open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)


def extract_features(input_path: str, output_bed: str, fmt: str = "auto", features=None):
    selected = set(features or ["gene", "transcript", "exon", "intron"])
    exons_by_tx = defaultdict(list)
    exon_idx_by_tx = defaultdict(int)
    rows = []

    for _line_no, record, error in iter_annotation(input_path, fmt):
        if error or record is None:
            continue
        kind = "transcript" if record.feature in {"mRNA", "lnc_RNA"} else record.feature
        if kind not in {"gene", "transcript", "exon"}:
            continue

        gid = gene_id(record)
        tid = transcript_id(record)
        gtype = gene_type(record)
        ttype = transcript_type(record)
        if kind == "gene" and "gene" in selected:
            rows.append(to_bed_row(record, "gene", gid, gid, "NA", gtype, "NA"))
        elif kind == "transcript" and "transcript" in selected:
            rows.append(to_bed_row(record, "transcript", tid, gid, tid, gtype, ttype))
        elif kind == "exon":
            key = (record.chrom, record.strand, gid, tid, gtype, ttype)
            exons_by_tx[key].append((record.start, record.end))
            exon_idx_by_tx[key] += 1
            exon_no = record.attrs.get("exon_number", str(exon_idx_by_tx[key]))
            if "exon" in selected:
                rows.append(to_bed_row(record, "exon", f"{tid}:exon{exon_no}", gid, tid, gtype, ttype))

    if "intron" in selected:
        for key, exons in exons_by_tx.items():
            chrom, strand, gid, tid, gtype, ttype = key
            exons_sorted = sorted(exons, key=lambda item: (item[0], item[1]))
            intron_idx = 0
            for left, right in zip(exons_sorted, exons_sorted[1:]):
                intron_start = left[1] + 1
                intron_end = right[0] - 1
                if intron_start <= intron_end:
                    intron_idx += 1
                    fake = AnnotationRecord(chrom, ".", "intron", intron_start, intron_end, ".", strand, ".", "", {})
                    rows.append(to_bed_row(fake, "intron", f"{tid}:intron{intron_idx}", gid, tid, gtype, ttype))

    rows.sort(key=lambda row: (row[0], int(row[1]), int(row[2]), row[3], row[4]))
    with open(output_bed, "w", encoding="utf-8") as out:
        for row in rows:
            out.write("\t".join(row) + "\n")


def extract_transcripts(input_path: str, output_bed: str, fmt: str = "auto"):
    extract_features(input_path, output_bed, fmt=fmt, features=["transcript"])

