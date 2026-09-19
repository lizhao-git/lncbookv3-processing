#!/usr/bin/env python3
import argparse
import csv


def esc_sql(s: str) -> str:
    return s.replace("'", "''")


def rgb_for_labels(labels: str) -> str:
    labels_lower = labels.lower()
    if "pathogenic" in labels_lower:
        return "255,0,0"
    if "benign" in labels_lower:
        return "0,128,0"
    if "risk factor" in labels_lower:
        return "255,140,0"
    if "drug response" in labels_lower:
        return "0,0,255"
    if "protective" in labels_lower:
        return "0,128,255"
    if "affects" in labels_lower:
        return "128,0,128"
    return "80,80,80"


def main():
    parser = argparse.ArgumentParser(
        description="Export annotation TSV to SQL and UCSC Genome Browser track BED"
    )
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
            "CREATE TABLE IF NOT EXISTS gtf_clinvar_annotation ("
            "chrom TEXT, start INTEGER, end INTEGER, feature_type TEXT, feature_id TEXT, "
            "gene_id TEXT, transcript_id TEXT, gene_type TEXT, transcript_type TEXT, strand TEXT, variant_id TEXT, variant_name TEXT, dbsnp_id TEXT, clinical_labels TEXT, "
            "allele_id TEXT, disease_name TEXT);\n"
        )
        out_bed.write(
            'track name="ClinVar_GTF_Site_Annotation" '
            'description="ClinVar labels on GTF gene/transcript/exon/intron" '
            'visibility=2 itemRgb="On"\n'
        )
        for r in reader:
            out_sql.write(
                "INSERT INTO gtf_clinvar_annotation "
                "(chrom, start, end, feature_type, feature_id, gene_id, transcript_id, gene_type, transcript_type, strand, variant_id, variant_name, dbsnp_id, clinical_labels, allele_id, disease_name) "
                "VALUES "
                f"('{esc_sql(r['chrom'])}', {int(r['start'])}, {int(r['end'])}, "
                f"'{esc_sql(r['feature_type'])}', '{esc_sql(r['feature_id'])}', "
                f"'{esc_sql(r['gene_id'])}', '{esc_sql(r['transcript_id'])}', "
                f"'{esc_sql(r.get('gene_type', 'NA'))}', '{esc_sql(r.get('transcript_type', 'NA'))}', "
                f"'{esc_sql(r['strand'])}', '{esc_sql(r['variant_id'])}', '{esc_sql(r.get('variant_name', 'NA'))}', '{esc_sql(r.get('dbsnp_id', 'NA'))}', '{esc_sql(r['clinical_labels'])}', "
                f"'{esc_sql(r.get('allele_id', 'NA'))}', '{esc_sql(r.get('disease_name', 'NA'))}');\n"
            )
            chrom = r["chrom"]
            start0 = int(r["start"]) - 1
            end1 = int(r["end"])
            name = f"{r['feature_type']}|{r['variant_id']}|{r['clinical_labels']}"
            strand = r["strand"] if r["strand"] in {"+", "-"} else "."
            rgb = rgb_for_labels(r["clinical_labels"])
            out_bed.write(
                f"{chrom}\t{start0}\t{end1}\t{name}\t0\t{strand}\t{start0}\t{end1}\t{rgb}\n"
            )
        out_sql.write("COMMIT;\n")


if __name__ == "__main__":
    main()
