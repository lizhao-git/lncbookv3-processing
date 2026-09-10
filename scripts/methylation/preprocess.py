#!/usr/bin/env python3
"""Stage A1-A3: normalize raw methylation data into a canonical 0-based,
4-column BED (``chrom  start  end  value``), with optional hg19->hg38 liftover.

Fixes legacy issues:
  * value column is unified to the 4th column (index 3) for every format;
  * BED coordinates are consistently 0-based half-open;
  * all paths are command-line parameters (no hard-coded ``/disk1/...``).

External tools required (on PATH):
  * ``bigWigToWig``, ``wig2bed``  (format=bigwig)
  * ``bigBedToBed``               (format=bigbed)
  * ``liftOver``                  (only when --chain is given)
"""

import argparse
import os
import shutil
import subprocess
import sys


def _run(cmd):
    print("+ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def _require(tool):
    path = shutil.which(tool)
    if path is None:
        sys.exit(f"error: required tool not found on PATH: {tool}")
    return path


def _pick_value(parts, preferred_index):
    """Return the methylation value as a float, falling back to the last
    parseable numeric column when the preferred index is missing/empty."""
    if preferred_index < len(parts):
        try:
            return float(parts[preferred_index])
        except ValueError:
            pass
    for tok in reversed(parts):
        try:
            return float(tok)
        except ValueError:
            continue
    return None


def _rewrite_bed(src, dst, value_index):
    """Rewrite a BED into the canonical 4-column form, extracting the value
    from ``value_index`` (0-based)."""
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            value = _pick_value(parts, value_index)
            if value is None:
                value = 0.0
            fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{value}\n")


def bigwig_to_bed(src, dst, value_index):
    _require("bigWigToWig")
    _require("wig2bed")
    os.makedirs(dst, exist_ok=True)
    for name in sorted(os.listdir(src)):
        stem, ext = os.path.splitext(name)
        if ext.lower() not in (".bw", ".bigwig"):
            continue
        wig = os.path.join(dst, stem + ".wig")
        _run(["bigWigToWig", os.path.join(src, name), wig])
        raw_bed = os.path.join(dst, stem + ".raw.bed")
        with open(raw_bed, "w") as out:
            subprocess.run(["wig2bed", "--zero-indexed"],
                           stdin=open(wig), stdout=out, check=True)
        _rewrite_bed(raw_bed, os.path.join(dst, stem + ".bed"), value_index)
        os.remove(raw_bed)
        os.remove(wig)


def bigbed_to_bed(src, dst, value_index):
    _require("bigBedToBed")
    os.makedirs(dst, exist_ok=True)
    for name in sorted(os.listdir(src)):
        stem, ext = os.path.splitext(name)
        if ext.lower() not in (".bb", ".bigbed"):
            continue
        raw_bed = os.path.join(dst, stem + ".raw.bed")
        _run(["bigBedToBed", os.path.join(src, name), raw_bed])
        _rewrite_bed(raw_bed, os.path.join(dst, stem + ".bed"), value_index)
        os.remove(raw_bed)


def bismark_to_bed(src, dst):
    """Parse Bismark CpG report files (``chrom pos strand meth unmeth ...``)
    into single-base BED with methylation rate = meth / (meth + unmeth)."""
    os.makedirs(dst, exist_ok=True)
    for name in sorted(os.listdir(src)):
        stem = name
        for suffix in ("_CpG_report.txt", ".CpG_report.txt", ".txt"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        bed = os.path.join(dst, stem + ".bed")
        with open(os.path.join(src, name)) as fin, open(bed, "w") as fout:
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos = parts[0], parts[1]
                try:
                    meth, unmeth = int(parts[3]), int(parts[4])
                    pos = int(pos)
                except ValueError:
                    continue
                total = meth + unmeth
                rate = (meth / total) if total else 0.0
                # 1-based CpG position -> 0-based half-open [pos-1, pos)
                fout.write(f"{chrom}\t{pos - 1}\t{pos}\t{rate}\n")


def bed_normalize(src, dst, value_index):
    os.makedirs(dst, exist_ok=True)
    for name in sorted(os.listdir(src)):
        if not name.endswith(".bed"):
            continue
        _rewrite_bed(os.path.join(src, name),
                     os.path.join(dst, name), value_index)


def liftover(src, dst, chain):
    lift = _require("liftOver")
    os.makedirs(dst, exist_ok=True)
    for name in sorted(os.listdir(src)):
        if not name.endswith(".bed"):
            continue
        stem = name[: -len(".bed")]
        _run([lift,
              os.path.join(src, name), chain,
              os.path.join(dst, stem + ".output.bed"),
              os.path.join(dst, stem + ".unmapped.bed")])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--format", required=True,
                    choices=["bigwig", "bigbed", "bismark", "bed"])
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--value-index", type=int, default=None,
                    help="0-based value column (default: 4 for bigwig, 3 otherwise)")
    ap.add_argument("--chain",
                    help="chain file for hg19->hg38 liftover (e.g. hg19ToHg38.over.chain.gz)")
    args = ap.parse_args()

    value_index = args.value_index if args.value_index is not None else \
        (4 if args.format == "bigwig" else 3)

    stage = os.path.join(args.output_dir, "bed")
    if args.format == "bigwig":
        bigwig_to_bed(args.input_dir, stage, value_index)
    elif args.format == "bigbed":
        bigbed_to_bed(args.input_dir, stage, value_index)
    elif args.format == "bismark":
        bismark_to_bed(args.input_dir, stage)
    else:  # bed
        bed_normalize(args.input_dir, stage, value_index)

    final_dir = stage
    if args.chain:
        hg38 = os.path.join(args.output_dir, "hg38_output")
        liftover(stage, hg38, args.chain)
        final_dir = hg38

    print(f"done: canonical BED dir = {final_dir}")


if __name__ == "__main__":
    main()
