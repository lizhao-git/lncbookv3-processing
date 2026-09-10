#!/usr/bin/env python3
import argparse


def main():
    parser = argparse.ArgumentParser(description="Format COSMIC bedtools intersections into annotation table")
    parser.add_argument("--input-intersections", required=True)
    parser.add_argument("--output-tsv", required=True)
    args = parser.parse_args()

    with open(args.input_intersections, "r", encoding="utf-8") as inp, open(args.output_tsv, "w", encoding="utf-8") as out:
        out.write(
            "chrom\tstart\tend\tfeature_type\tfeature_id\tgene_id\ttranscript_id\tgene_type\ttranscript_type\tstrand\tvariant_id\tclinical_labels\tcosmic_id\tdisease_name\tfathmm_mkl_score\n"
        )
        for raw in inp:
            line = raw.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 18:
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
            labels = cols[14]
            cosmic_id = cols[15]
            disease_name = cols[16]
            fathmm_mkl_score = cols[17]
            out.write(
                f"{chrom}\t{start_1}\t{end_1}\t{feature_type}\t{feature_id}\t{gene_id}\t{tx_id}\t{gene_type}\t{transcript_type}\t{strand}\t{variant_id}\t{labels}\t{cosmic_id}\t{disease_name}\t{fathmm_mkl_score}\n"
            )


if __name__ == "__main__":
    main()
