#!/usr/bin/env python3
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description="Retain SmProt mappings entirely and uniquely within lncRNA transcripts")
    parser.add_argument("--input-intersections", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    contained_rows = []
    transcripts_by_protein = defaultdict(set)
    total_rows = 0

    with open(args.input_intersections, "r", encoding="utf-8") as inp:
        for raw in inp:
            line = raw.rstrip("\n")
            if not line:
                continue
            total_rows += 1
            cols = line.split("\t")
            if len(cols) < 16:
                continue

            tx_chrom = cols[0]
            tx_start = int(cols[1])
            tx_end = int(cols[2])
            feature_type = cols[3]
            feature_id = cols[4]
            gene_id = cols[5]
            transcript_id = cols[6]
            strand = cols[7]
            gene_type = cols[8] if len(cols) > 8 else "NA"
            transcript_type = cols[9] if len(cols) > 9 else "NA"

            protein_chrom = cols[10]
            protein_start = int(cols[11])
            protein_end = int(cols[12])
            protein_id = cols[13]
            evidence = cols[14] if len(cols) > 14 else "NA"
            source = cols[15] if len(cols) > 15 else "SmProt"

            if feature_type in {"gene", "transcript"}:
                continue
            if tx_chrom != protein_chrom:
                continue
            if not (protein_start >= tx_start and protein_end <= tx_end):
                continue

            key = (protein_id, transcript_id)
            transcripts_by_protein[protein_id].add(transcript_id)
            contained_rows.append((
                key,
                tx_chrom,
                protein_start,
                protein_end,
                protein_id,
                evidence,
                source,
                transcript_id,
                gene_id,
                gene_type,
                transcript_type,
                strand,
                tx_start,
                tx_end,
                feature_id,
            ))

    unique_proteins = {pid for pid, txs in transcripts_by_protein.items() if len(txs) == 1}
    retained = [row for row in contained_rows if row[0][0] in unique_proteins]

    retained.sort(key=lambda r: (r[1], r[2], r[3], r[4]))

    with open(args.output_tsv, "w", encoding="utf-8") as out:
        out.write(
            "chrom\tstart\tend\tprotein_id\tevidence\tsource\ttranscript_id\tgene_id\tgene_type\ttranscript_type\tstrand\ttranscript_start\ttranscript_end\tmapping_type\n"
        )
        for row in retained:
            _key, chrom, start0, end1, protein_id, evidence, source, tx_id, gene_id, gene_type, transcript_type, strand, tx_start0, tx_end1, _feature_id = row
            out.write(
                f"{chrom}\t{start0}\t{end1}\t{protein_id}\t{evidence}\t{source}\t{tx_id}\t{gene_id}\t{gene_type}\t{transcript_type}\t{strand}\t{tx_start0}\t{tx_end1}\tentirely_within_unique_transcript\n"
            )

    retained_proteins = {r[4] for r in retained}
    retained_genes = {r[8] for r in retained}
    retained_transcripts = {r[7] for r in retained}

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"total_overlap_rows\t{total_rows}\n")
        rep.write(f"contained_overlap_rows\t{len(contained_rows)}\n")
        rep.write(f"retained_rows\t{len(retained)}\n")
        rep.write(f"retained_small_proteins\t{len(retained_proteins)}\n")
        rep.write(f"retained_lncRNA_genes\t{len(retained_genes)}\n")
        rep.write(f"retained_lncRNA_transcripts\t{len(retained_transcripts)}\n")


if __name__ == "__main__":
    main()
