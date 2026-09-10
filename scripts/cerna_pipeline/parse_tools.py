#!/usr/bin/env python3
"""Parse the three target-prediction tools' outputs into a single normalized
interaction format::

    miRNA  target  score  energy  start  end

Input conventions (matching the LncBook 2022 legacy files):

* miRanda   : ``>miRNA  target  score  energy  start  end``
* TargetScan: ``target  miRNA  start  end``
* RNAhybrid : ``miRNA  target  start  end  energy  pvalue``
              (single file, or a directory of ``*_new_out.txt`` split files)
"""

import argparse
import os


def parse_miranda(src, dst):
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            miRNA = parts[0].lstrip(">")
            target = parts[1]
            score, energy = parts[2], parts[3]
            start, end = parts[4], parts[5]
            fout.write(f"{miRNA}\t{target}\t{score}\t{energy}\t{start}\t{end}\n")


def parse_targetscan(src, dst):
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            target, miRNA = parts[0], parts[1]
            start, end = parts[2], parts[3]
            fout.write(f"{miRNA}\t{target}\t-\t-\t{start}\t{end}\n")


def _rnahybrid_paths(src):
    if os.path.isdir(src):
        return [os.path.join(src, n) for n in sorted(os.listdir(src))
                if n.endswith(".txt")]
    return [src]


def parse_rnahybrid(src, dst):
    with open(dst, "w") as fout:
        for path in _rnahybrid_paths(src):
            with open(path) as fin:
                for line in fin:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    miRNA, target = parts[0], parts[1]
                    start, end = parts[2], parts[3]
                    energy = parts[4]
                    fout.write(f"{miRNA}\t{target}\t-\t{energy}\t{start}\t{end}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--miranda")
    ap.add_argument("--targetscan")
    ap.add_argument("--rnahybrid")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    if args.miranda:
        parse_miranda(args.miranda, os.path.join(args.output_dir, "miranda.tsv"))
    if args.targetscan:
        parse_targetscan(args.targetscan, os.path.join(args.output_dir, "targetscan.tsv"))
    if args.rnahybrid:
        parse_rnahybrid(args.rnahybrid, os.path.join(args.output_dir, "rnahybrid.tsv"))
    print(f"done: normalized tables -> {args.output_dir}")


if __name__ == "__main__":
    main()
