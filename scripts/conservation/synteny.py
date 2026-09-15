#!/usr/bin/env python3
import argparse
import os
import sys

if __package__ is None or __package__ == "":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from conservation.common import read_tsv, write_tsv


def key_for(row, prefix):
    return (
        row.get(f"{prefix}_gene_id") or row.get(f"{prefix}_gene") or row.get(f"{prefix}gene_id"),
        row.get("target_species") or row.get("species"),
    )


def call_synteny(homology_rows, anchors):
    anchor_set = set()
    for row in anchors:
        human = row.get("human_gene_id") or row.get("human_gene")
        species = row.get("target_species") or row.get("species")
        target = row.get("target_gene_id") or row.get("target_gene")
        if human and species and target:
            anchor_set.add((human, species, target))

    out = []
    for row in homology_rows:
        human = row.get("human_gene_id")
        species = row.get("target_species")
        target = row.get("target_gene_id")
        rec = dict(row)
        rec["syntenic_support"] = "1" if (human, species, target) in anchor_set else "0"
        out.append(rec)
    return out


def main():
    ap = argparse.ArgumentParser(description="Annotate homology rows with optional synteny anchor support")
    ap.add_argument("--homology", required=True)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    write_tsv(args.output, call_synteny(read_tsv(args.homology), read_tsv(args.anchors)))


if __name__ == "__main__":
    main()
