#!/usr/bin/env python3
import argparse
import csv


def esc_sql(s: str) -> str:
    return s.replace("'", "''")


def main():
    parser = argparse.ArgumentParser(description="Export COSMIC annotation TSV to SQL and UCSC Genome Browser track BED")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-sql", required=True)
    parser.add_argument("--output-ucsc-bed", required=True)
    args = parser.parse_args()

    # Single streaming pass: both outputs are written row by row so memory
    # stays constant regardless of the annotation table size.
    with open(args.input_tsv, "r", encoding="utf-8") as fh, open(
        args.output_sql, "w", encoding="utf-8"
    ) as out_sql, open(args.output_ucsc_bed, "w", encoding="utf-8") as out_bed:
        reader = csv.DictReader(fh, delimiter="\t")
        out_sql.write("BEGIN TRANSACTION;\n")
        out_sql.write(
            "CREATE TABLE IF NOT EXISTS gtf_cosmic_annotation ("
            "chrom TEXT, start INTEGER, end INTEGER, feature_type TEXT, feature_id TEXT, gene_id TEXT, transcript_id TEXT, gene_type TEXT, transcript_type TEXT, strand TEXT, variant_id TEXT, clinical_labels TEXT, cosmic_id TEXT, disease_name TEXT, fathmm_mkl_score REAL);\n"
        )
        out_bed.write(
            'track name="COSMIC_GTF_Site_Annotation" description="COSMIC pathogenic variants on GTF gene/transcript/exon/intron" visibility=2 itemRgb="On"\n'
        )
        for r in reader:
            out_sql.write(
                "INSERT INTO gtf_cosmic_annotation (chrom, start, end, feature_type, feature_id, gene_id, transcript_id, gene_type, transcript_type, strand, variant_id, clinical_labels, cosmic_id, disease_name, fathmm_mkl_score) VALUES "
                f"('{esc_sql(r['chrom'])}', {int(r['start'])}, {int(r['end'])}, '{esc_sql(r['feature_type'])}', '{esc_sql(r['feature_id'])}', '{esc_sql(r['gene_id'])}', '{esc_sql(r['transcript_id'])}', '{esc_sql(r.get('gene_type', 'NA'))}', '{esc_sql(r.get('transcript_type', 'NA'))}', '{esc_sql(r['strand'])}', '{esc_sql(r['variant_id'])}', '{esc_sql(r['clinical_labels'])}', '{esc_sql(r['cosmic_id'])}', '{esc_sql(r['disease_name'])}', {float(r['fathmm_mkl_score'])});\n"
            )
            chrom = r['chrom']
            start0 = int(r['start']) - 1
            end1 = int(r['end'])
            name = f"{r['feature_type']}|{r['cosmic_id']}|{r['disease_name']}"
            out_bed.write(f"{chrom}\t{start0}\t{end1}\t{name}\t0\t{r['strand']}\t{start0}\t{end1}\t255,0,0\n")
        out_sql.write("COMMIT;\n")


if __name__ == '__main__':
    main()
