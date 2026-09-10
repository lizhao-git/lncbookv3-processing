#!/usr/bin/env python3
import argparse


def main():
    parser = argparse.ArgumentParser(description="Format GWAS Catalog bedtools intersections into annotation table")
    parser.add_argument("--input-intersections", required=True)
    parser.add_argument("--output-tsv", required=True)
    args = parser.parse_args()

    with open(args.input_intersections, "r", encoding="utf-8") as inp, open(args.output_tsv, "w", encoding="utf-8") as out:
        out.write(
            "chrom\tstart\tend\tfeature_type\tfeature_id\tgene_id\ttranscript_id\tgene_type\ttranscript_type\tstrand\tvariant_id\tassociation_label\tgwas_variant_id\tallele_id\ttrait_name\tefo_trait\tefo_id\tp_value\trisk_allele_frequency\tor_or_beta\tci_95_text\n"
        )
        for raw in inp:
            line = raw.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 24:
                continue
            chrom = cols[0]
            start_1 = int(cols[1]) + 1
            end_1 = int(cols[2])
            feature_type = cols[3]
            feature_id = cols[4]
            gene_id = cols[5]
            tx_id = cols[6]
            strand = cols[7]
            gene_type = cols[8]
            transcript_type = cols[9]
            if feature_type in {"gene", "transcript"}:
                continue
            variant_id = cols[13]
            association_label = cols[14]
            gwas_variant_id = cols[15]
            trait_name = cols[16]
            p_value = cols[17]
            allele_id = cols[18]
            efo_trait = cols[19]
            efo_id = cols[20]
            risk_allele_frequency = cols[21]
            or_or_beta = cols[22]
            ci_95_text = cols[23]
            out.write(
                f"{chrom}\t{start_1}\t{end_1}\t{feature_type}\t{feature_id}\t{gene_id}\t{tx_id}\t{gene_type}\t{transcript_type}\t{strand}\t{variant_id}\t{association_label}\t{gwas_variant_id}\t{allele_id}\t{trait_name}\t{efo_trait}\t{efo_id}\t{p_value}\t{risk_allele_frequency}\t{or_or_beta}\t{ci_95_text}\n"
            )


if __name__ == "__main__":
    main()
