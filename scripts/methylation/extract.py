#!/usr/bin/env python3
"""Stage A4-A5: per-gene mean methylation over the gene body / promoter
regions, then combine samples into gene x sample matrices.

Fixes legacy issues:
  * single-pass, path-parameterized, Python 3;
  * correct 0-based overlap (no off-by-one triple-condition);
  * unified value column (canonical 4-column BED);
  * no 10k-gene chunking needed (removes the misleading ``*_0_10000`` names);
  * chromosome intervals are loaded once per sample (legacy re-read the
    chromosome file for every gene).

Input:  a directory of per-sample 4-column BED files (from ``preprocess.py``).
Output: ``{prefix}_body`` and ``{prefix}_promoter`` matrices with two header
        lines (sample names, then sample classes), plus per-sample value files.
"""

import argparse
import json
import os
import re
import sys

from common import load_gene_info, mean_over_region


def read_bed_by_chrom(path):
    """Load a 4-column BED into ``{chrom: sorted[(start, end, value), ...]}``."""
    by_chrom = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start, end, value = int(parts[1]), int(parts[2]), float(parts[3])
            except ValueError:
                continue
            by_chrom.setdefault(chrom, []).append((start, end, value))
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda t: t[0])
    return by_chrom


def extract_sample(bed_path, genes):
    """Return ``{gene_id: mean_value}`` for the genes overlapping the BED."""
    by_chrom = read_bed_by_chrom(bed_path)
    result = {}
    for gid, (chrom, s0, e0) in genes.items():
        intervals = by_chrom.get(chrom)
        if intervals is None:
            continue
        mean = mean_over_region(intervals, s0, e0)
        if mean is not None:
            result[gid] = mean
    return result


def classify_sample(sample, class_map, regex_map):
    if class_map and sample in class_map:
        return class_map[sample]
    for cls, rx in (regex_map or {}).items():
        if rx.search(sample):
            return cls
    return "NA"


def parse_class_args(class_map_file, class_regex):
    class_map = None
    if class_map_file:
        with open(class_map_file) as fh:
            class_map = json.load(fh)
    regex_map = None
    if class_regex:
        regex_map = {}
        for spec in class_regex.split(";"):
            if "=" not in spec:
                continue
            cls, pattern = spec.split("=", 1)
            regex_map[cls.strip()] = re.compile(pattern.strip())
    return class_map, regex_map


def write_matrix(path, gene_ids, per_sample, sample_names, classes):
    with open(path, "w") as fh:
        fh.write("Gene\t" + "\t".join(sample_names) + "\n")
        fh.write("Gene\t" + "\t".join(classes) + "\n")
        for gid in gene_ids:
            vals = []
            for name in sample_names:
                v = per_sample[name].get(gid)
                vals.append(str(v) if v is not None else "0.0")
            fh.write(gid + "\t" + "\t".join(vals) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed-dir", required=True,
                    help="directory of per-sample 4-column BED files")
    ap.add_argument("--gene-info", required=True,
                    help="gene body reference (gene_id chrom start end, 1-based closed)")
    ap.add_argument("--promoter-info",
                    help="promoter reference (same format); omit to skip promoter")
    ap.add_argument("--prefix", required=True,
                    help="output prefix, e.g. GSE109875")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--class-map",
                    help="JSON file mapping sample stem -> class label")
    ap.add_argument("--class-regex",
                    help="class regexes, e.g. 'case=MB|tumor;control=Normal'")
    args = ap.parse_args()

    genes = load_gene_info(args.gene_info)
    promoters = load_gene_info(args.promoter_info) if args.promoter_info else None
    class_map, regex_map = parse_class_args(args.class_map, args.class_regex)

    sample_files = sorted(f for f in os.listdir(args.bed_dir)
                          if f.endswith(".bed") and not f.endswith(".unmapped.bed"))
    if not sample_files:
        sys.exit(f"error: no .bed files under {args.bed_dir}")

    gene_ids = sorted(genes)
    os.makedirs(args.output_dir, exist_ok=True)

    body_per_sample = {}
    promoter_per_sample = {}
    sample_names = []
    classes = []

    for fname in sample_files:
        stem = fname[: -len(".bed")]
        bed_path = os.path.join(args.bed_dir, fname)
        print(f"extracting {stem}", file=sys.stderr)

        body_values = extract_sample(bed_path, genes)
        body_per_sample[stem] = body_values
        with open(os.path.join(args.output_dir, f"{stem}_body.txt"), "w") as fh:
            for gid in gene_ids:
                if gid in body_values:
                    fh.write(f"{gid}\t{body_values[gid]}\n")

        if promoters is not None:
            promoter_values = extract_sample(bed_path, promoters)
            promoter_per_sample[stem] = promoter_values
            with open(os.path.join(args.output_dir, f"{stem}_promoter.txt"), "w") as fh:
                for gid in gene_ids:
                    if gid in promoter_values:
                        fh.write(f"{gid}\t{promoter_values[gid]}\n")

        sample_names.append(stem)
        classes.append(classify_sample(stem, class_map, regex_map))

    write_matrix(os.path.join(args.output_dir, f"{args.prefix}_body"),
                 gene_ids, body_per_sample, sample_names, classes)
    if promoters is not None:
        write_matrix(os.path.join(args.output_dir, f"{args.prefix}_promoter"),
                     gene_ids, promoter_per_sample, sample_names, classes)

    print(f"done: {len(sample_files)} samples -> {args.output_dir}")


if __name__ == "__main__":
    main()
