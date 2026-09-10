#!/usr/bin/env python3
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description="Summarize retained SmProt-lncRNA mappings")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-summary", required=True)
    args = parser.parse_args()

    row_count = 0
    proteins = set()
    genes = set()
    transcripts = set()

    with open(args.input_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row_count += 1
            proteins.add(row.get("protein_id", ""))
            genes.add(row.get("gene_id", ""))
            transcripts.add(row.get("transcript_id", ""))

    proteins.discard("")
    genes.discard("")
    transcripts.discard("")

    with open(args.output_summary, "w", encoding="utf-8") as out:
        out.write("metric\tvalue\n")
        out.write(f"retained_rows\t{row_count}\n")
        out.write(f"retained_small_proteins\t{len(proteins)}\n")
        out.write(f"retained_lncRNA_genes\t{len(genes)}\n")
        out.write(f"retained_lncRNA_transcripts\t{len(transcripts)}\n")


if __name__ == "__main__":
    main()
