import heapq
import os
import re
import tempfile
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
    if lower.endswith(".gz"):
        lower = lower[:-3]
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


EXPLICIT_UTR_FEATURES = {"UTR", "five_prime_utr", "three_prime_utr", "5UTR", "3UTR"}


def _resolve_utr_side(strand: str, start: int, end: int, cds_min, cds_max, explicit=None) -> str:
    """Classify a UTR segment as 'utr5' / 'utr3' / 'utr' (unclassifiable).

    ``start``/``end`` are 1-based GTF coordinates. ``explicit`` may already state
    the prime end (GENCODE-style ``five_prime_utr`` / ``three_prime_utr``); when
    it is unknown the segment is positioned relative to the CDS span: on the
    + strand lower coordinates are 5', on the - strand lower coordinates are 3'.
    """
    if explicit == "five_prime_utr":
        return "utr5"
    if explicit == "three_prime_utr":
        return "utr3"
    if cds_min is None or cds_max is None:
        return "utr"
    if end < cds_min:
        before_cds = True
    elif start > cds_max:
        before_cds = False
    else:
        return "utr"
    if strand == "-":
        return "utr3" if before_cds else "utr5"
    return "utr5" if before_cds else "utr3"


_FEATURE_CHUNK_ROWS = 250_000


def _feature_row_key(row):
    return (row[0], int(row[1]), int(row[2]), row[3], row[4])


def _flush_feature_chunk(rows, chunk_paths, work_dir):
    """Sort one buffered chunk of feature rows and spill it to a temp file."""
    rows.sort(key=_feature_row_key)
    fd, path = tempfile.mkstemp(prefix=".features_chunk_", suffix=".tsv", dir=work_dir)
    with os.fdopen(fd, "w", encoding="utf-8") as chunk_fh:
        for row in rows:
            chunk_fh.write("\t".join(row) + "\n")
    chunk_paths.append(path)


def extract_features(input_path: str, output_bed: str, fmt: str = "auto", features=None):
    """Extract gene/transcript/exon/intron/CDS/UTR intervals into a BED-like table.

    ``features`` accepts the lowercase tokens ``gene``, ``transcript``, ``exon``,
    ``intron``, ``cds`` and ``utr`` (``utr`` expands to both ``utr5`` and
    ``utr3``); ``utr5`` / ``utr3`` select a single UTR end. UTR intervals are
    taken from explicit GENCODE-style ``UTR`` / ``five_prime_utr`` /
    ``three_prime_utr`` rows when present, and otherwise derived from the
    exon minus CDS layout. Transcripts without a CDS have no UTR.
    """
    selected = set(features or ["gene", "transcript", "exon", "intron"])
    want_utr = bool({"utr", "utr5", "utr3"} & selected)
    want_utr5 = want_utr and ("utr" in selected or "utr5" in selected)
    want_utr3 = want_utr and ("utr" in selected or "utr3" in selected)

    exons_by_tx = defaultdict(list)
    cds_by_tx = defaultdict(list)
    explicit_utr_by_tx = defaultdict(list)
    exon_idx_by_tx = defaultdict(int)
    tx_meta = {}
    # Output rows are sorted by (chrom, start, end, kind, id). Large
    # annotations (e.g. full GENCODE) yield millions of rows, far beyond the
    # process memory budget, so the buffer is flushed in sorted chunks and the
    # chunks are merged at the end instead of holding every row in memory.
    chunk_paths = []
    work_dir = os.path.dirname(os.path.abspath(output_bed)) or "."
    rows = []

    def emit(row):
        rows.append(row)
        if len(rows) >= _FEATURE_CHUNK_ROWS:
            _flush_feature_chunk(rows, chunk_paths, work_dir)
            rows.clear()

    for _line_no, record, error in iter_annotation(input_path, fmt):
        if error or record is None:
            continue
        feature = record.feature
        kind = "transcript" if feature in {"mRNA", "lnc_RNA"} else feature
        gid = gene_id(record)
        tid = transcript_id(record)
        gtype = gene_type(record)
        ttype = transcript_type(record)
        key = (record.chrom, record.strand, gid, tid)
        tx_meta[key] = (gid, gtype, ttype)

        if kind == "gene":
            if "gene" in selected:
                emit(to_bed_row(record, "gene", gid, gid, "NA", gtype, "NA"))
        elif kind == "transcript":
            if "transcript" in selected:
                emit(to_bed_row(record, "transcript", tid, gid, tid, gtype, ttype))
        elif kind == "exon":
            exons_by_tx[key].append((record.start, record.end))
            exon_idx_by_tx[key] += 1
            exon_no = record.attrs.get("exon_number", str(exon_idx_by_tx[key]))
            if "exon" in selected:
                emit(to_bed_row(record, "exon", f"{tid}:exon{exon_no}", gid, tid, gtype, ttype))
        elif feature == "CDS":
            cds_by_tx[key].append((record.start, record.end))
            if "cds" in selected:
                emit(to_bed_row(record, "CDS", f"{tid}:CDS", gid, tid, gtype, ttype))
        elif feature in EXPLICIT_UTR_FEATURES:
            explicit_utr_by_tx[key].append((record.start, record.end, feature))

    if "intron" in selected:
        for key, exons in exons_by_tx.items():
            chrom, strand, gid, tid = key
            _, gtype, ttype = tx_meta.get(key, (gid, "NA", "NA"))
            exons_sorted = sorted(exons, key=lambda item: (item[0], item[1]))
            intron_idx = 0
            for left, right in zip(exons_sorted, exons_sorted[1:]):
                intron_start = left[1] + 1
                intron_end = right[0] - 1
                if intron_start <= intron_end:
                    intron_idx += 1
                    fake = AnnotationRecord(chrom, ".", "intron", intron_start, intron_end, ".", strand, ".", "", {})
                    emit(to_bed_row(fake, "intron", f"{tid}:intron{intron_idx}", gid, tid, gtype, ttype))

    if want_utr:
        for key, exons in exons_by_tx.items():
            chrom, strand, gid, tid = key
            _, gtype, ttype = tx_meta.get(key, (gid, "NA", "NA"))
            cds_list = cds_by_tx.get(key, [])
            if cds_list:
                cds_min = min(s for s, _ in cds_list)
                cds_max = max(e for _, e in cds_list)
            else:
                cds_min = cds_max = None

            segments = []  # (start, end, explicit_prime)
            explicit = explicit_utr_by_tx.get(key, [])
            if explicit:
                for (s, e, fname) in explicit:
                    prime = fname if fname in {"five_prime_utr", "three_prime_utr"} else None
                    if fname == "5UTR":
                        prime = "five_prime_utr"
                    elif fname == "3UTR":
                        prime = "three_prime_utr"
                    segments.append((s, e, prime))
            elif cds_min is not None:
                for (e_start, e_end) in sorted(exons):
                    if e_end < cds_min:
                        segments.append((e_start, e_end, None))
                    elif e_start > cds_max:
                        segments.append((e_start, e_end, None))
                    else:
                        if e_start < cds_min:
                            segments.append((e_start, cds_min - 1, None))
                        if e_end > cds_max:
                            segments.append((cds_max + 1, e_end, None))

            for (s, e, prime) in segments:
                if s > e:
                    continue
                which = _resolve_utr_side(strand, s, e, cds_min, cds_max, explicit=prime)
                if which == "utr5" and not want_utr5:
                    continue
                if which == "utr3" and not want_utr3:
                    continue
                label = {"utr5": "5UTR", "utr3": "3UTR"}.get(which, "UTR")
                fake = AnnotationRecord(chrom, ".", label, s, e, ".", strand, ".", "", {})
                emit(to_bed_row(fake, label, f"{tid}:{label}", gid, tid, gtype, ttype))

    with open(output_bed, "w", encoding="utf-8") as out:
        if not chunk_paths:
            rows.sort(key=_feature_row_key)
            for row in rows:
                out.write("\t".join(row) + "\n")
            return
        if rows:
            _flush_feature_chunk(rows, chunk_paths, work_dir)
            rows.clear()
        try:
            readers = [open(path, "r", encoding="utf-8") for path in chunk_paths]
            try:
                merged = heapq.merge(
                    *readers, key=lambda line: _feature_row_key(line.split("\t"))
                )
                for line in merged:
                    out.write(line)
            finally:
                for reader in readers:
                    reader.close()
        finally:
            for path in chunk_paths:
                try:
                    os.remove(path)
                except OSError:
                    pass


def extract_transcripts(input_path: str, output_bed: str, fmt: str = "auto"):
    extract_features(input_path, output_bed, fmt=fmt, features=["transcript"])

