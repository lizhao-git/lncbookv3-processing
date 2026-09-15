#!/usr/bin/env python3
"""Post-process pslMap conservation results.

This is the Python replacement for the legacy group_process_filter_lnc.R and
group_process_filter_pc.R scripts. It consumes paired query/target FASTA,
mapped BED12/group files and transcript context tables, then emits one
per-species statistics table.
"""

import argparse
import os
import sys
from collections import defaultdict

_SCRIPTS_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if __package__ is None or __package__ == "":
    sys.path.insert(0, _SCRIPTS_ROOT)

# BED12 block/interval handling lives in the shared format_convert layer.
_SHARED_LIB_ROOT = os.path.join(_SCRIPTS_ROOT, "format_convert")
if _SHARED_LIB_ROOT not in sys.path:
    sys.path.insert(0, _SHARED_LIB_ROOT)

from genomic_intervals.bed12 import BED12_FIELDS, BlockRecord, interval_union_length, overlap_fraction_of_smaller
from conservation.common import parse_int, read_fasta, read_tsv, write_tsv


def base_id(segment_id):
    parts = segment_id.split("_")
    if len(parts) <= 1:
        return segment_id
    return "_".join(parts[:-1])


def matches(query, target):
    return sum(1 for q, t in zip(query.upper(), target.upper()) if q == t)


def read_block_table(path):
    rows = read_tsv(path, header=False, names=BED12_FIELDS)
    return {row["ids"]: BlockRecord.from_row(row) for row in rows}


def choose_duplicate_drop(segment_ids, data_by_segment, length_by_segment, overlap_cutoff):
    drops = set()
    grouped = defaultdict(list)
    for segment_id in segment_ids:
        grouped[base_id(segment_id)].append(segment_id)

    for _gid, members in grouped.items():
        members = sorted(members)
        for i, left_id in enumerate(members):
            if left_id in drops:
                continue
            left = length_by_segment.get(left_id)
            if left is None:
                continue
            for right_id in members[i + 1:]:
                if right_id in drops:
                    continue
                right = length_by_segment.get(right_id)
                if right is None or not left.span_intersects(right):
                    continue
                frac = overlap_fraction_of_smaller(left.intervals(), right.intervals())
                if frac < overlap_cutoff:
                    continue
                left_score = parse_int(data_by_segment[left_id]["matched_length"])
                right_score = parse_int(data_by_segment[right_id]["matched_length"])
                if left_score >= right_score:
                    drops.add(right_id)
                else:
                    drops.add(left_id)
                    break
    return drops


def nearby_groups(segment_ids, group_by_segment, max_distance):
    grouped = defaultdict(list)
    for segment_id in segment_ids:
        grouped[base_id(segment_id)].append(segment_id)

    renamed = {}
    modified = set()
    for gid, members in grouped.items():
        members = sorted(members)
        if len(members) == 1:
            renamed[members[0]] = gid
            continue

        labels = list(range(len(members)))
        for i, left_id in enumerate(members):
            left = group_by_segment.get(left_id)
            if left is None:
                continue
            for j in range(i + 1, len(members)):
                right_id = members[j]
                right = group_by_segment.get(right_id)
                if right is None:
                    continue
                if left.strand == right.strand and left.distance_within(right, max_distance):
                    keep = min(labels[i], labels[j])
                    labels[i] = keep
                    labels[j] = keep
        label_map = {old: idx + 1 for idx, old in enumerate(sorted(set(labels)))}
        for segment_id, label in zip(members, labels):
            grouped_id = f"{gid}_{label_map[label]}"
            renamed[segment_id] = grouped_id
            modified.add(grouped_id)
    return renamed, modified


def covered_exons(length_record, exon_record):
    if length_record is None or exon_record is None:
        return [], 0
    q_intervals = length_record.intervals()
    covered = []
    for idx, exon_interval in enumerate(exon_record.intervals(), start=1):
        if any(max(qs, exon_interval[0]) < min(qe, exon_interval[1]) for qs, qe in q_intervals):
            covered.append(idx)
    return covered, len(exon_record.intervals())


def aggregate_group(grouped_id, members, rows, length_by_segment, exon_by_id):
    representative = max(members, key=lambda sid: parse_int(rows[sid].get("matched_length")))
    out = dict(rows[representative])
    out["ids"] = base_id(representative) if "_" in representative else representative
    out["grouped_ids"] = grouped_id
    out["mapped_length"] = sum(parse_int(rows[sid].get("mapped_length")) for sid in members)
    out["matched_length"] = sum(parse_int(rows[sid].get("matched_length")) for sid in members)
    out["match_ratio"] = safe_ratio(out["matched_length"], out["mapped_length"])

    intervals = []
    covered_union = set()
    exon_count = 0
    for sid in members:
        length_record = length_by_segment.get(sid)
        if length_record:
            intervals.extend(length_record.intervals())
        covered, count = covered_exons(length_record, exon_by_id.get(base_id(sid)))
        exon_count = max(exon_count, count)
        covered_union.update(covered)
    out["true_mapped_length"] = interval_union_length(intervals)
    out["exon_counts"] = exon_count
    out["exon_covered"] = ",".join(str(i) for i in sorted(covered_union))
    out["exon_covered_num"] = len(covered_union)
    out["exon_covered_ratio"] = safe_ratio(len(covered_union), exon_count)
    return out


def safe_ratio(num, den):
    den = float(den or 0)
    if den <= 0:
        return ""
    return float(num or 0) / den


def merge_by_key(rows, table, key, prefix=""):
    lookup = {row.get(key): row for row in table if row.get(key)}
    out = []
    for row in rows:
        merged = dict(row)
        hit = lookup.get(row.get(key))
        if hit:
            for hit_key, value in hit.items():
                if hit_key == key:
                    continue
                out_key = f"{prefix}{hit_key}" if prefix else hit_key
                if out_key in merged:
                    out_key = f"{prefix or 'context_'}{hit_key}"
                merged[out_key] = value
        out.append(merged)
    return out


def run_postprocess(args):
    query = read_fasta(args.query_fasta)
    target = read_fasta(args.target_fasta)
    ids = sorted(query)
    target_ids = sorted(target)
    if ids != target_ids:
        raise SystemExit("query/target FASTA IDs do not match after sorting")

    data_by_segment = {}
    for sid in ids:
        qseq = query[sid]
        tseq = target[sid]
        if len(qseq) != len(tseq):
            raise SystemExit(f"mapped lengths differ for {sid}: {len(qseq)} vs {len(tseq)}")
        matched = matches(qseq, tseq)
        mapped = len(qseq)
        data_by_segment[sid] = {
            "segment_id": sid,
            "ids": base_id(sid),
            "mapped_length": mapped,
            "matched_length": matched,
            "match_ratio": safe_ratio(matched, mapped),
        }

    group_by_segment = read_block_table(args.group_info)
    length_by_segment = read_block_table(args.lengths)
    exon_by_id = read_block_table(args.exon_positions)

    drops = choose_duplicate_drop(ids, data_by_segment, length_by_segment, args.overlap_cutoff)
    kept = [sid for sid in ids if sid not in drops]
    renamed, modified = nearby_groups(kept, group_by_segment, args.nearby_distance)

    groups = defaultdict(list)
    for sid in kept:
        groups[renamed[sid]].append(sid)

    rows = [aggregate_group(gid, members, data_by_segment, length_by_segment, exon_by_id)
            for gid, members in sorted(groups.items())]

    if args.target_info:
        target_info = read_tsv(args.target_info, header=True)
        rows = merge_by_key(rows, target_info, "ids")
    if args.target_addition:
        addition = read_tsv(args.target_addition, header=True)
        key = "target_transcript_id" if args.mode == "lnc" else "target_transcript_id"
        rows = merge_by_key(rows, addition, key)
    if args.transcript_lengths:
        lengths = read_tsv(args.transcript_lengths, header=False, names=["ids", "transcript_length"])
        rows = merge_by_key(rows, lengths, "ids", "")
    if args.transcript_meta:
        meta_names = [
            "chromosome", "ids", "strand", "transcript_type", "transcript_name",
            "transcript_support_level", "geneID", "gene_name", "gene_type",
            "gene_classification", "transcript_classification",
        ]
        meta = read_tsv(args.transcript_meta, header=False, names=meta_names)
        rows = merge_by_key(rows, meta, "ids", "")

    for row in rows:
        row["modified_group"] = "1" if row.get("grouped_ids") in modified else "0"
        row["dropped_duplicate_segments"] = len(drops)
        row["transcript_coverage"] = safe_ratio(
            parse_int(row.get("true_mapped_length")),
            parse_int(row.get("transcript_length")),
        )

    if args.mode == "pc":
        rows = add_pc_metrics(rows, query, target, args)

    preferred = [
        "ids", "grouped_ids", "geneID", "gene_name", "transcript_length",
        "true_mapped_length", "mapped_length", "matched_length", "match_ratio",
        "transcript_coverage", "exon_counts", "exon_covered", "exon_covered_num",
        "exon_covered_ratio", "modified_group", "dropped_duplicate_segments",
    ]
    fields = preferred + sorted({key for row in rows for key in row if key not in preferred})
    write_tsv(args.output, rows, fields)
    write_tsv(args.report, [{
        "mode": args.mode,
        "segments_input": len(ids),
        "segments_dropped": len(drops),
        "segments_kept": len(kept),
        "groups_output": len(rows),
    }])


def add_pc_metrics(rows, query, target, args):
    if not args.cds_positions:
        return rows
    cds_rows = read_tsv(args.cds_positions, header=False, names=["ids", "CDS_start", "CDS_end"])
    cds_by_id = {row["ids"]: row for row in cds_rows}
    for row in rows:
        sid = row.get("segment_id") or row.get("ids")
        cds = cds_by_id.get(sid) or cds_by_id.get(row.get("ids"))
        if not cds:
            continue
        start = parse_int(cds.get("CDS_start"), 1) - 1
        end = parse_int(cds.get("CDS_end"))
        qseq = query.get(sid, "")
        tseq = target.get(sid, "")
        if end <= start or not qseq or not tseq:
            continue
        q_sub = qseq[start:end]
        t_sub = tseq[start:end]
        cds_matched = matches(q_sub, t_sub)
        cds_mapped = min(len(q_sub), len(t_sub))
        row["CDS_mapped_length"] = cds_mapped
        row["CDS_matched_length"] = cds_matched
        row["CDS_match_ratio"] = safe_ratio(cds_matched, cds_mapped)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["lnc", "pc"], required=True)
    ap.add_argument("--query-fasta", required=True)
    ap.add_argument("--target-fasta", required=True)
    ap.add_argument("--group-info", required=True)
    ap.add_argument("--lengths", required=True)
    ap.add_argument("--exon-positions", required=True)
    ap.add_argument("--target-info")
    ap.add_argument("--target-addition")
    ap.add_argument("--transcript-lengths")
    ap.add_argument("--transcript-meta")
    ap.add_argument("--cds-positions")
    ap.add_argument("--overlap-cutoff", type=float, default=0.5)
    ap.add_argument("--nearby-distance", type=int, default=100000)
    ap.add_argument("--output", required=True)
    ap.add_argument("--report", required=True)
    run_postprocess(ap.parse_args())


if __name__ == "__main__":
    main()
