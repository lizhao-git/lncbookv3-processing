#!/usr/bin/env python3
import argparse
import os
import sys
from collections import defaultdict

if __package__ is None or __package__ == "":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from conservation.common import parse_number, read_tsv, write_tsv


def score(row):
    return (
        parse_number(row.get("true_mapped_length")),
        parse_number(row.get("match_ratio")),
        parse_number(row.get("transcript_coverage")),
        parse_number(row.get("exon_covered_ratio")),
    )


def gene_id(row):
    return row.get("geneID") or row.get("gene_id") or row.get("ids") or ""


def target_gene_id(row):
    for key in ("target_gene_id", "target_geneID", "target_gene", "gene_id_target", "target_extra_gene_id"):
        if row.get(key):
            return row[key]
    return row.get("target_transcript_id") or row.get("target_ids") or ""


def build_homology(rows):
    evidence = []
    for row in rows:
        if row.get("is_conserved") not in ("1", "true", "True", "yes", "YES"):
            continue
        record = {
            "human_gene_id": gene_id(row),
            "human_transcript_id": row.get("ids", ""),
            "target_species": row.get("target_species", ""),
            "target_gene_id": target_gene_id(row),
            "target_transcript_id": row.get("target_transcript_id") or row.get("target_ids") or "",
            "true_mapped_length": row.get("true_mapped_length", ""),
            "match_ratio": row.get("match_ratio", ""),
            "transcript_coverage": row.get("transcript_coverage", ""),
            "exon_covered_ratio": row.get("exon_covered_ratio", ""),
            "evidence_group": row.get("grouped_ids", ""),
        }
        evidence.append(record)

    groups = defaultdict(list)
    for row in evidence:
        key = (row["human_gene_id"], row["target_species"], row["target_gene_id"])
        groups[key].append(row)

    best = []
    for (human_gene, species, target_gene), members in sorted(groups.items()):
        chosen = max(members, key=score)
        out = dict(chosen)
        out["human_gene_id"] = human_gene
        out["target_species"] = species
        out["target_gene_id"] = target_gene
        out["supporting_transcript_count"] = len({m["human_transcript_id"] for m in members})
        out["supporting_segment_count"] = len(members)
        best.append(out)
    return evidence, best


def main():
    ap = argparse.ArgumentParser(description="Collapse transcript-level conservation to gene homology")
    ap.add_argument("--conservation", required=True)
    ap.add_argument("--output-all", required=True)
    ap.add_argument("--output-best", required=True)
    args = ap.parse_args()

    evidence, best = build_homology(read_tsv(args.conservation))
    fields = [
        "human_gene_id", "human_transcript_id", "target_species",
        "target_gene_id", "target_transcript_id", "true_mapped_length",
        "match_ratio", "transcript_coverage", "exon_covered_ratio",
        "evidence_group",
    ]
    write_tsv(args.output_all, evidence, fields)
    write_tsv(args.output_best, best)


if __name__ == "__main__":
    main()
