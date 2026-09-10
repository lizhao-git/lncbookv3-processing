#!/usr/bin/env python3
import argparse
import csv


def esc_sql(value: str) -> str:
    return value.replace("'", "''")


def main():
    parser = argparse.ArgumentParser(description="Export GWAS Catalog annotation TSV to SQL and UCSC Genome Browser track BED")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-sql", required=True)
    parser.add_argument("--output-ucsc-bed", required=True)
    args = parser.parse_args()

    with open(args.output_sql, "w", encoding="utf-8") as out_sql:
        out_sql.write("BEGIN TRANSACTION;\n")
        out_sql.write(
            "CREATE TABLE IF NOT EXISTS gtf_gwas_catalog_annotation ("
            "chrom TEXT, start INTEGER, end INTEGER, feature_type TEXT, feature_id TEXT, gene_id TEXT, transcript_id TEXT, gene_type TEXT, transcript_type TEXT, strand TEXT, variant_id TEXT, association_label TEXT, gwas_variant_id TEXT, allele_id TEXT, trait_name TEXT, efo_trait TEXT, efo_id TEXT, p_value REAL, risk_allele_frequency TEXT, or_or_beta TEXT, ci_95_text TEXT);\n"
        )
        with open(args.input_tsv, "r", encoding="utf-8") as in_fh:
            reader = csv.DictReader(in_fh, delimiter="\t")
            for row in reader:
                out_sql.write(
                    "INSERT INTO gtf_gwas_catalog_annotation (chrom, start, end, feature_type, feature_id, gene_id, transcript_id, gene_type, transcript_type, strand, variant_id, association_label, gwas_variant_id, allele_id, trait_name, efo_trait, efo_id, p_value, risk_allele_frequency, or_or_beta, ci_95_text) VALUES "
                    f"('{esc_sql(row['chrom'])}', {int(row['start'])}, {int(row['end'])}, '{esc_sql(row['feature_type'])}', '{esc_sql(row['feature_id'])}', '{esc_sql(row['gene_id'])}', '{esc_sql(row['transcript_id'])}', '{esc_sql(row.get('gene_type', 'NA'))}', '{esc_sql(row.get('transcript_type', 'NA'))}', '{esc_sql(row['strand'])}', '{esc_sql(row['variant_id'])}', '{esc_sql(row['association_label'])}', '{esc_sql(row['gwas_variant_id'])}', '{esc_sql(row.get('allele_id', 'NA'))}', '{esc_sql(row['trait_name'])}', '{esc_sql(row.get('efo_trait', 'NA'))}', '{esc_sql(row.get('efo_id', 'NA'))}', {float(row['p_value'])}, '{esc_sql(row.get('risk_allele_frequency', 'NA'))}', '{esc_sql(row.get('or_or_beta', 'NA'))}', '{esc_sql(row.get('ci_95_text', 'NA'))}');\n"
                )
        out_sql.write("COMMIT;\n")

    with open(args.output_ucsc_bed, "w", encoding="utf-8") as out_bed:
        out_bed.write(
            'track name="GWAS_Catalog_GTF_Site_Annotation" description="GWAS Catalog genome-wide significant variants on GTF gene/transcript/exon/intron" visibility=2 itemRgb="On"\n'
        )
        with open(args.input_tsv, "r", encoding="utf-8") as in_fh:
            reader = csv.DictReader(in_fh, delimiter="\t")
            for row in reader:
                chrom = row["chrom"]
                start0 = int(row["start"]) - 1
                end1 = int(row["end"])
                name = f"{row['feature_type']}|{row['gwas_variant_id']}|{row['trait_name']}"
                out_bed.write(f"{chrom}\t{start0}\t{end1}\t{name}\t0\t{row['strand']}\t{start0}\t{end1}\t0,102,204\n")


if __name__ == "__main__":
    main()
