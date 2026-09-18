#!/usr/bin/env python3
import argparse
import csv
from collections import Counter


def main():
    parser = argparse.ArgumentParser(description="Summarize feature-level ClinVar annotations")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-summary", required=True)
    args = parser.parse_args()

    LABELS = ["Pathogenic", "Benign", "Affects", "Drug response", "Protective", "Risk factor"]
    # Most specific location class first, so the summary reads top-down.
    PREFERRED = ["CDS", "5UTR", "3UTR", "exon", "intron", "gene", "transcript"]

    by_feature = Counter()
    by_feature_label = Counter()

    with open(args.input_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            feature = row.get("feature_type", "")
            labels_raw = row.get("clinical_labels", "")
            labels = labels_raw.split("|") if labels_raw else []
            by_feature[feature] += 1
            for lb in labels:
                by_feature_label[(feature, lb)] += 1

    ordered = [f for f in PREFERRED if f in by_feature]
    ordered += sorted(f for f in by_feature if f not in PREFERRED)

    with open(args.output_summary, "w", encoding="utf-8") as out:
        out.write("feature_type\tlabel\tcount\n")
        for feature in ordered:
            out.write(f"{feature}\tALL\t{by_feature.get(feature, 0)}\n")
            for lb in LABELS:
                out.write(f"{feature}\t{lb}\t{by_feature_label.get((feature, lb), 0)}\n")


if __name__ == "__main__":
    main()
