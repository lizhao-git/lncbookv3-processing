#!/usr/bin/env python3
"""Experimental validation, final interaction table, and disease annotation.

Reimplements the legacy ``code/step6..step8`` downstream steps:

* experimental ceRNA mapping (LncRNAWiki + HGNC -> gene-miRNA pairs)
* final interaction table (transcript -> gene, with a header row)
* HMDD disease annotation (miRNA-disease -> lncRNA-disease)

Reference file formats (LncBook 2022 legacy):
  * transcript_id_alias.txt : ``gene_id  ENST...;NR_...;XR_...``
  * HGNC_symbol_transcript.txt : ``symbol  alias(|-sep)  refseq``
  * transcript_gene.txt     : ``transcript  gene``
  * LncRNAWiki2.0-ceRNA.txt : ``lncRNA_symbol  miRNA;miRNA;...``
  * HMDD alldata_v2.txt     : ``...  miRNA  disease``
"""

import argparse
import os
import re


def build_transcript_gene(path):
    """transcript_id_alias.txt -> {transcript_id(no version): gene_id}."""
    id_convert = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene = parts[0]
            for alias in parts[1].split(";"):
                if alias.startswith(("ENST", "NR_", "XR_")):
                    id_convert[alias.split(".")[0]] = gene
    return id_convert


def build_hgnc(path):
    """HGNC_symbol_transcript.txt -> {symbol: [transcript ids]}."""
    hgnc = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            symbol = parts[0]
            content = []
            if parts[1] != "!!!":
                content.extend(a for a in parts[1].split("|") if a.startswith("ENST"))
            if parts[2] != "!!!":
                content.append(parts[2])
            hgnc[symbol] = content
    return hgnc


def experiment(lncrnawiki, hgnc, id_convert, out_path):
    """LncRNAWiki + HGNC -> gene-miRNA experimental pairs."""
    all_ceRNA = {}
    with open(lncrnawiki) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            symbol = parts[0]
            if symbol not in hgnc:
                continue
            for transcript in hgnc[symbol]:
                if transcript not in id_convert:
                    continue
                gene = id_convert[transcript]
                for mir in parts[1].split(";"):
                    if mir.startswith("miR"):
                        mir = "hsa-" + mir
                    if mir.startswith("hsa"):
                        all_ceRNA.setdefault(gene, []).append(mir)
    with open(out_path, "w") as fw:
        for gene in all_ceRNA:
            for mir in sorted(set(all_ceRNA[gene])):
                fw.write(f"{gene}\t{mir}\n")


def final_table(three_overlap, transcript_gene, out_path):
    """three_overlap -> final interaction table with a header."""
    tran_gene = {}
    with open(transcript_gene) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                tran_gene[parts[0]] = parts[1]
    with open(three_overlap) as fin, open(out_path, "w") as fw:
        fw.write("GeneID\tLncID\tmiRNA\tScore\tEnergy\tStart\tEnd\n")
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            lnc = parts[0]
            gene = tran_gene.get(lnc, lnc)
            fw.write(gene + "\t" + "\t".join(parts) + "\n")


def _strip_arm(mir):
    for arm in ("-3p", "-5p"):
        if arm in mir:
            return mir.split(arm)[0]
    return mir


def disease_annotate(hmdd, final_interaction, out_path):
    """HMDD -> lncRNA-disease annotation."""
    miRNA = {}
    with open(hmdd, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            miRNA.setdefault(parts[1], []).append(parts[2])
    miRNA_used = {m: d for m, d in miRNA.items() if len(d) >= 2}

    lnc = {}
    with open(final_interaction) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] == "GeneID":
                continue
            mir = _strip_arm(parts[2])
            if mir in miRNA_used:
                lnc.setdefault(parts[1], []).append(mir)

    lnc_disease = {}
    for l, mirs in lnc.items():
        if len(mirs) < 2:
            continue
        for m in mirs:
            if m in miRNA_used:
                lnc_disease.setdefault(l, []).extend(miRNA_used[m])

    with open(out_path, "w") as fw:
        for l in lnc_disease:
            fw.write(f"{l}\t{';'.join(lnc_disease[l])}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lncrnawiki", required=True)
    ap.add_argument("--hgnc", required=True)
    ap.add_argument("--transcript-id-alias", required=True)
    ap.add_argument("--transcript-gene", required=True)
    ap.add_argument("--hmdd", required=True)
    ap.add_argument("--three-overlap", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    id_convert = build_transcript_gene(args.transcript_id_alias)
    hgnc = build_hgnc(args.hgnc)

    experiment(args.lncrnawiki, hgnc, id_convert,
               os.path.join(args.output_dir, "ceRNA_experiment_out.txt"))

    final_path = os.path.join(args.output_dir, "LncBook_v2_interaction_out.txt")
    final_table(args.three_overlap, args.transcript_gene, final_path)

    disease_annotate(args.hmdd, final_path,
                     os.path.join(args.output_dir, "LncBook_v2_interaction_out_HMDD.txt"))

    print(f"done: experimental + final + disease annotation -> {args.output_dir}")


if __name__ == "__main__":
    main()
