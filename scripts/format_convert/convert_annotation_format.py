#!/usr/bin/env python3
"""Convert annotation files between GTF and GFF3.

GTF -> GFF3 keeps every original attribute (including gene_id/transcript_id)
and adds ID/Parent according to the GFF3 hierarchy: gene features get
ID=gene_id, transcript features get ID=transcript_id plus Parent=gene_id,
and child features (exon, CDS, ...) get Parent=transcript_id so they hang
off their transcript rather than the gene. Feature types, coordinates,
score, strand and frame/phase are passed through unchanged.

GFF3 -> GTF restores gene_id/transcript_id by classifying each feature from
its Parent chain (transcripts are direct children of gene-like features,
everything else is a child of its transcript) and keeps all other
attributes. GFF3 directive comments (``##...``) are dropped when writing
GTF, and the ``##gff-version 3`` line is added when writing GFF3.
"""
import argparse

from genomic_intervals.annotation import (
    detect_annotation_format,
    gff3_escape,
    gff3_unescape,
    parse_gff3_attrs,
    parse_gtf_attrs,
)
from genomic_intervals.io import open_text

TRANSCRIPT_FEATURES = {"transcript", "mRNA", "lnc_RNA", "ncRNA", "tRNA", "rRNA"}
GENE_FEATURES = {"gene", "pseudogene"}


def is_gene_like(feature: str) -> bool:
    return feature in GENE_FEATURES or feature.endswith("_gene")


def convert_gtf_to_gff3(input_path: str, output_path: str) -> dict:
    stats = {
        "input_records": 0,
        "output_records": 0,
        "gene_features": 0,
        "transcript_features": 0,
        "child_features": 0,
        "records_missing_gene_id": 0,
        "records_missing_transcript_id": 0,
    }
    with open_text(input_path) as (fh, _compression), open(output_path, "w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("##gff-version"):
                    continue
                out.write(line + "\n")
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                raise SystemExit(f"{input_path}: line {line_no}: expected 9 columns, got {len(cols)}")
            chrom, source, feature, start, end, score, strand, frame, attrs_text = cols
            stats["input_records"] += 1
            attrs = parse_gtf_attrs(attrs_text)
            gid = attrs.get("gene_id")
            tid = attrs.get("transcript_id")

            pairs = []
            if is_gene_like(feature):
                role = "gene"
                if gid:
                    pairs.append(("ID", gid))
                else:
                    stats["records_missing_gene_id"] += 1
            elif feature in TRANSCRIPT_FEATURES:
                role = "transcript"
                if tid:
                    pairs.append(("ID", tid))
                    if gid:
                        pairs.append(("Parent", gid))
                    else:
                        stats["records_missing_gene_id"] += 1
                else:
                    stats["records_missing_transcript_id"] += 1
                    if gid:
                        pairs.append(("ID", gid))
                    else:
                        stats["records_missing_gene_id"] += 1
            else:
                role = "child"
                if tid:
                    pairs.append(("Parent", tid))
                    if not gid:
                        stats["records_missing_gene_id"] += 1
                elif gid:
                    pairs.append(("Parent", gid))
                    stats["records_missing_transcript_id"] += 1
                else:
                    stats["records_missing_gene_id"] += 1

            stats[f"{role}_features"] += 1

            all_pairs = pairs + [(key, value) for key, value in attrs.items()]
            encoded = ";".join(f"{key}={gff3_escape(value)}" for key, value in all_pairs)
            out.write("\t".join([chrom, source, feature, start, end, score, strand, frame, encoded]) + "\n")
            stats["output_records"] += 1
    return stats


def build_feature_maps(input_path: str):
    """Map feature ID to its (first) Parent ID and to its feature type."""
    parent_of = {}
    type_of = {}
    with open_text(input_path) as (fh, _compression):
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                continue
            attrs = parse_gff3_attrs(cols[8])
            id_ = attrs.get("ID")
            if not id_:
                continue
            type_of[id_] = cols[2]
            parent = attrs.get("Parent")
            if parent:
                parent_of[id_] = parent.split(",")[0]
    return parent_of, type_of


def resolve_ids(parent_id: str, parent_of: dict, type_of: dict):
    """Walk up the Parent chain collecting ``(gene_id, transcript_id)``."""
    gene = None
    transcript = None
    seen = set()
    current = parent_id
    while current and current not in seen:
        seen.add(current)
        if is_gene_like(type_of.get(current, "")):
            if gene is None:
                gene = current
        elif transcript is None:
            transcript = current
        current = parent_of.get(current)
    return gene, transcript


def escape_gtf_value(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def convert_gff3_to_gtf(input_path: str, output_path: str) -> dict:
    stats = {
        "input_records": 0,
        "output_records": 0,
        "gene_features": 0,
        "transcript_features": 0,
        "child_features": 0,
        "records_missing_gene_id": 0,
        "records_missing_transcript_id": 0,
    }
    parent_of, type_of = build_feature_maps(input_path)
    with open_text(input_path) as (fh, _compression), open(output_path, "w", encoding="utf-8") as out:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("##"):
                    continue
                out.write(line + "\n")
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                raise SystemExit(f"{input_path}: line {line_no}: expected 9 columns, got {len(cols)}")
            chrom, source, feature, start, end, score, strand, frame, attrs_text = cols
            stats["input_records"] += 1
            attrs = parse_gff3_attrs(attrs_text)
            id_ = attrs.get("ID")
            parent = attrs.get("Parent")
            if parent:
                parent = parent.split(",")[0]

            gid = None
            tid = None
            if parent:
                gid, ancestor_tid = resolve_ids(parent, parent_of, type_of)
                if is_gene_like(type_of.get(parent, "")):
                    role = "transcript"
                    tid = id_ or ancestor_tid
                else:
                    role = "child"
                    tid = ancestor_tid or id_
            elif feature in TRANSCRIPT_FEATURES:
                role = "transcript"
                tid = id_
            else:
                role = "gene"
                gid = id_
            stats[f"{role}_features"] += 1
            if not gid:
                stats["records_missing_gene_id"] += 1
            if not tid and role != "gene":
                stats["records_missing_transcript_id"] += 1

            pairs = [("gene_id", gid or ".")]
            if role != "gene":
                pairs.append(("transcript_id", tid or "."))
            for key, value in attrs.items():
                if key in ("ID", "Parent", "gene_id", "transcript_id"):
                    continue
                pairs.append((key, gff3_unescape(value)))

            encoded = " ".join(f'{key} "{escape_gtf_value(value)}";' for key, value in pairs)
            out.write("\t".join([chrom, source, feature, start, end, score, strand, frame, encoded]) + "\n")
            stats["output_records"] += 1
    return stats


def write_report(path: str, direction: str, stats: dict):
    with open(path, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"direction\t{direction}\n")
        for key, value in stats.items():
            rep.write(f"{key}\t{value}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert GTF <-> GFF3 annotation files")
    parser.add_argument("--input-annotation", required=True)
    parser.add_argument("--output-annotation", required=True)
    parser.add_argument("--from", dest="from_format", choices=["auto", "gtf", "gff3"], default="auto")
    parser.add_argument("--to", dest="to_format", choices=["gtf", "gff3"], required=True)
    parser.add_argument("--report")
    args = parser.parse_args()

    source_format = detect_annotation_format(args.input_annotation, args.from_format)
    if source_format == args.to_format:
        raise SystemExit(f"Input is already {source_format}; nothing to convert.")

    if args.to_format == "gff3":
        stats = convert_gtf_to_gff3(args.input_annotation, args.output_annotation)
        direction = "gtf_to_gff3"
    else:
        stats = convert_gff3_to_gtf(args.input_annotation, args.output_annotation)
        direction = "gff3_to_gtf"

    if args.report:
        write_report(args.report, direction, stats)


if __name__ == "__main__":
    main()
