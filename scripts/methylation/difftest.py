#!/usr/bin/env python3
"""Stage B (previously missing): differential methylation testing.

Reads a gene x sample matrix produced by ``extract.py`` (two header lines:
sample names, then sample classes) and calls significantly differentially
methylated genes between case and control samples.

Methods:
  * ``wilcox``     : two-sided Mann-Whitney U + Benjamini-Hochberg FDR.
                     A gene is significant when p < --pval (default 0.02).
                     Direction = High/Low by median_case vs median_control,
                     and the winning side must exceed --min-median (0.2).
  * ``foldchange`` : TCGA/liver-style consistency test. A gene is "High" when
                     every case value > every control value, fold change
                     (min_case / max_control) >= --fc and the winning side
                     exceeds --min-abs; symmetric for "Low".

Outputs (per region):
  * ``{prefix}_{region}_test.txt``       full wilcox results
        gene  p  median_case  median_control
  * ``{prefix}_{region}_sig.txt``        significant calls (both methods)
        gene  direction  evidence
"""

import argparse
import os
import sys

from common import bh_fdr, mann_whitney_u, median


def read_matrix(path):
    with open(path) as fh:
        lines = fh.readlines()
    if len(lines) < 3:
        sys.exit(f"error: matrix has too few lines: {path}")
    samples = lines[0].rstrip("\n").split("\t")[1:]
    classes = lines[1].rstrip("\n").split("\t")[1:]
    rows = {}
    for line in lines[2:]:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        try:
            rows[parts[0]] = [float(v) for v in parts[1:]]
        except ValueError:
            continue
    return samples, classes, rows


def split_indices(classes, case, control):
    case_idx = [i for i, c in enumerate(classes) if c == case]
    ctrl_idx = [i for i, c in enumerate(classes) if c == control]
    if not case_idx or not ctrl_idx:
        sys.exit(f"error: need >=1 sample for both '{case}' and '{control}' "
                 f"(found classes: {sorted(set(classes))})")
    return case_idx, ctrl_idx


def wilcox(rows, case_idx, ctrl_idx, pval, min_median):
    genes, pvalues, med_cases, med_ctrls = [], [], [], []
    for gid, vals in rows.items():
        a = [vals[i] for i in case_idx]
        b = [vals[i] for i in ctrl_idx]
        _, p = mann_whitney_u(a, b)
        genes.append(gid)
        pvalues.append(p)
        med_cases.append(median(a))
        med_ctrls.append(median(b))

    qvalues = bh_fdr(pvalues)
    results = []
    for gid, p, q, mc, mn in zip(genes, pvalues, qvalues, med_cases, med_ctrls):
        direction = None
        if mc > mn and mc > min_median:
            direction = "High"
        elif mn > mc and mn > min_median:
            direction = "Low"
        evidence = (f"p={p:.6g};fdr={q:.6g};"
                    f"median_case={mc:.6g};median_control={mn:.6g}")
        results.append((gid, p, q, mc, mn, direction, evidence))
    return results


def foldchange(rows, case_idx, ctrl_idx, fc, min_abs):
    results = []
    for gid, vals in rows.items():
        a = [vals[i] for i in case_idx]
        b = [vals[i] for i in ctrl_idx]
        min_a, max_a = min(a), max(a)
        min_b, max_b = min(b), max(b)
        direction = None
        evidence = ""
        if min_a > max_b:  # all cases higher than all controls
            ratio = min_a / (max_b + 1e-9)
            if ratio >= fc and max_a > min_abs:
                direction = "High"
                evidence = (f"FC={ratio:.6g};min_case={min_a:.6g};"
                            f"max_control={max_b:.6g}")
        elif max_a < min_b:  # all cases lower than all controls
            ratio = max_b / (min_a + 1e-9)
            if ratio >= fc and max_b > min_abs:
                direction = "Low"
                evidence = (f"FC={ratio:.6g};max_case={max_a:.6g};"
                            f"min_control={min_b:.6g}")
        results.append((gid, direction, evidence))
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--region", choices=["body", "promoter"], required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--method", choices=["wilcox", "foldchange"], default="wilcox")
    ap.add_argument("--case-label", default="case")
    ap.add_argument("--control-label", default="control")
    ap.add_argument("--pval", type=float, default=0.02)
    ap.add_argument("--min-median", type=float, default=0.2,
                    help="wilcox: winning side median must exceed this")
    ap.add_argument("--fc", type=float, default=2.0,
                    help="foldchange: minimum fold change")
    ap.add_argument("--min-abs", type=float, default=20.0,
                    help="foldchange: winning side must exceed this")
    args = ap.parse_args()

    samples, classes, rows = read_matrix(args.matrix)
    case_idx, ctrl_idx = split_indices(classes, args.case_label, args.control_label)
    os.makedirs(args.output_dir, exist_ok=True)

    prefix = args.prefix
    region = args.region
    test_path = os.path.join(args.output_dir, f"{prefix}_{region}_test.txt")
    sig_path = os.path.join(args.output_dir, f"{prefix}_{region}_sig.txt")

    if args.method == "wilcox":
        results = wilcox(rows, case_idx, ctrl_idx, args.pval, args.min_median)
        with open(test_path, "w") as fh:
            for gid, p, q, mc, mn, _d, _e in results:
                fh.write(f"{gid}\t{p:.6g}\t{mc:.6g}\t{mn:.6g}\n")
        with open(sig_path, "w") as fh:
            for gid, p, q, mc, mn, direction, evidence in results:
                if direction is not None and p < args.pval:
                    fh.write(f"{gid}\t{direction}\t{evidence}\n")
    else:
        results = foldchange(rows, case_idx, ctrl_idx, args.fc, args.min_abs)
        with open(sig_path, "w") as fh:
            for gid, direction, evidence in results:
                if direction is not None:
                    fh.write(f"{gid}\t{direction}\t{evidence}\n")

    print(f"done: {sig_path}")


if __name__ == "__main__":
    main()
