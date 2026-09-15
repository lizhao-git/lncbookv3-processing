#!/usr/bin/env python3
import argparse
import os
import sys

if __package__ is None or __package__ == "":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from conservation.common import read_tsv, write_tsv


def load_species_order(path):
    rows = read_tsv(path)
    order = {}
    label = {}
    for idx, row in enumerate(rows):
        species = row.get("species") or row.get("target_species")
        if not species:
            continue
        rank = row.get("rank") or row.get("age_rank") or idx
        try:
            rank = int(rank)
        except ValueError:
            rank = idx
        order[species] = rank
        label[species] = row.get("clade") or row.get("age") or species
    return order, label


def assign_gene_age(homology_rows, species_order_path):
    order, label = load_species_order(species_order_path)
    by_gene = {}
    for row in homology_rows:
        gene = row.get("human_gene_id")
        species = row.get("target_species")
        if not gene or not species:
            continue
        rank = order.get(species, 10**9)
        current = by_gene.get(gene)
        if current is None or rank < current["rank"]:
            by_gene[gene] = {
                "human_gene_id": gene,
                "oldest_detected_species": species,
                "oldest_detected_clade": label.get(species, species),
                "rank": rank,
                "supporting_target_genes": set(),
                "supporting_transcripts": set(),
                "species": set(),
            }
        target = row.get("target_gene_id")
        tx = row.get("human_transcript_id")
        by_gene[gene]["species"].add(species)
        if target:
            by_gene[gene]["supporting_target_genes"].add(target)
        if tx:
            by_gene[gene]["supporting_transcripts"].add(tx)

    out = []
    for gene, rec in sorted(by_gene.items()):
        out.append({
            "human_gene_id": gene,
            "oldest_detected_species": rec["oldest_detected_species"],
            "oldest_detected_clade": rec["oldest_detected_clade"],
            "species_count": len(rec["species"]),
            "supporting_target_genes": ",".join(sorted(rec["supporting_target_genes"])),
            "supporting_transcripts": ",".join(sorted(rec["supporting_transcripts"])),
        })
    return out


def main():
    ap = argparse.ArgumentParser(description="Assign lncRNA gene age from the oldest homologous species")
    ap.add_argument("--homology", required=True)
    ap.add_argument("--species-order", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    rows = assign_gene_age(read_tsv(args.homology), args.species_order)
    write_tsv(args.output, rows)


if __name__ == "__main__":
    main()
