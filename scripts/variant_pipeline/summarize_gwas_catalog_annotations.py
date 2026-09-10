#!/usr/bin/env python3
import argparse
import csv
from collections import Counter

LABEL = "Genome-wide significant"


def main():
    parser = argparse.ArgumentParser(description="Summarize feature-level GWAS Catalog annotations")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-summary", required=True)
    args = parser.parse_args()

    by_feature = Counter()
    by_feature_label = Counter()

    with open(args.input_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            feature = row.get("feature_type", "")
            label = row.get("association_label", "")
            by_feature[feature] += 1
            if label:
                by_feature_label[(feature, label)] += 1

    with open(args.output_summary, "w", encoding="utf-8") as out:
        out.write("feature_type\tlabel\tcount\n")
        for feature in ["gene", "transcript", "exon", "intron"]:
            out.write(f"{feature}\tALL\t{by_feature.get(feature, 0)}\n")
            out.write(f"{feature}\t{LABEL}\t{by_feature_label.get((feature, LABEL), 0)}\n")


if __name__ == "__main__":
    main()
