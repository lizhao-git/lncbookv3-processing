#!/usr/bin/env python3
import argparse
import csv


def esc_sql(value: str) -> str:
    return value.replace("'", "''")


def main():
    parser = argparse.ArgumentParser(description="Export SmProt annotation TSV to SQL and UCSC BED")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-sql", required=True)
    parser.add_argument("--output-ucsc-bed", required=True)
    args = parser.parse_args()

    with open(args.output_sql, "w", encoding="utf-8") as out_sql:
        out_sql.write("BEGIN TRANSACTION;\n")
        out_sql.write(
            "CREATE TABLE IF NOT EXISTS gtf_smprot_annotation ("
            "chrom TEXT, start INTEGER, end INTEGER, protein_id TEXT, evidence TEXT, source TEXT, transcript_id TEXT, gene_id TEXT, gene_type TEXT, transcript_type TEXT, strand TEXT, transcript_start INTEGER, transcript_end INTEGER, mapping_type TEXT);\n"
        )

        with open(args.input_tsv, "r", encoding="utf-8") as in_fh:
            reader = csv.DictReader(in_fh, delimiter="\t")
            for row in reader:
                out_sql.write(
                    "INSERT INTO gtf_smprot_annotation (chrom, start, end, protein_id, evidence, source, transcript_id, gene_id, gene_type, transcript_type, strand, transcript_start, transcript_end, mapping_type) VALUES "
                    f"('{esc_sql(row['chrom'])}', {int(row['start'])}, {int(row['end'])}, '{esc_sql(row['protein_id'])}', '{esc_sql(row['evidence'])}', '{esc_sql(row['source'])}', '{esc_sql(row['transcript_id'])}', '{esc_sql(row['gene_id'])}', '{esc_sql(row.get('gene_type', 'NA'))}', '{esc_sql(row.get('transcript_type', 'NA'))}', '{esc_sql(row['strand'])}', {int(row['transcript_start'])}, {int(row['transcript_end'])}, '{esc_sql(row['mapping_type'])}');\n"
                )
        out_sql.write("COMMIT;\n")

    with open(args.output_ucsc_bed, "w", encoding="utf-8") as out_bed:
        out_bed.write(
            'track name="SmProt_lncRNA_Annotation" description="SmProt small proteins entirely and uniquely within lncRNA transcripts" visibility=2 itemRgb="On"\n'
        )
        with open(args.input_tsv, "r", encoding="utf-8") as in_fh:
            reader = csv.DictReader(in_fh, delimiter="\t")
            for row in reader:
                chrom = row["chrom"]
                start0 = int(row["start"])
                end1 = int(row["end"])
                name = f"{row['protein_id']}|{row['transcript_id']}|{row['gene_id']}"
                strand = row["strand"] if row["strand"] in {"+", "-"} else "."
                out_bed.write(f"{chrom}\t{start0}\t{end1}\t{name}\t0\t{strand}\t{start0}\t{end1}\t102,0,204\n")


if __name__ == "__main__":
    main()
