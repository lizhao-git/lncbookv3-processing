#!/usr/bin/env python3
"""Run miRanda / RNAhybrid / TargetScan on the prepared sequence inputs and
emit the three normalized prediction files consumed by the downstream ceRNA
pipeline (parse_tools -> predict_ceRNA -> experiment_annotate).

Outputs (under --output-dir):
  miranda/miRanda_interaction_out_select_v2.txt   >miRNA target score energy start end
  rnahybrid/RNAhybrid_out.txt                     miRNA target start end energy pvalue
  targetscan/TargetScan_combine_v2.txt            target miRNA start end

The three start/end columns are 1-based inclusive binding sites on the
lncRNA (target) sequence, so the downstream sliding-window refinement
compares positions in one consistent coordinate system.

miRanda and RNAhybrid are invoked as external binaries (see Dockerfile);
TargetScan follows the legacy simplified logic: the miRNA seed (2-8 nt,
7-mer) is reverse-complemented and matched directly against the lncRNA
sequence.
"""

import argparse
import os
import subprocess
import sys

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def run(cmd, out_path=None):
    print("+ " + " ".join(cmd), file=sys.stderr)
    if out_path:
        with open(out_path, "w") as fh:
            subprocess.run(cmd, stdout=fh, check=True)
    else:
        subprocess.run(cmd, check=True)


def parse_miranda(raw, out):
    """miRanda hit block header -> ``>miRNA target score energy start end``.

    Columns (tab-separated): 0=>miRNA 1=target 2=total_score 3=total_energy
    4=max_score 5=max_energy 6=qry_start 7=qry_end 8=sbj_start 9=sbj_end ...
    start/end are the target (lncRNA) binding site.
    """
    with open(raw) as fin, open(out, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            mirna, target, score, energy = p[0], p[1], p[2], p[3]
            start, end = p[8], p[9]
            fout.write(f"{mirna}\t{target}\t{score}\t{energy}\t{start}\t{end}\n")


def parse_rnahybrid(raw, out, pvalue_cutoff):
    """RNAhybrid human-readable block -> ``miRNA target start end energy pvalue``.

    Block layout (0-based line index from the ``target:`` line):
      0 target: X, 1 length: N, 2 miRNA : X, 3 length: n, 4 (blank),
      5 mfe: .., 6 p-value: .., 7 (blank), 8 position N, 9 target 5' ..,
      10 alignment line, 11 miRNA 3' ..
    """
    lines = open(raw, encoding="utf-8", errors="replace").read().splitlines()
    with open(out, "w") as fout:
        for i, line in enumerate(lines):
            if not line.startswith("target:"):
                continue
            try:
                target = line.split(":", 1)[1].strip()
                mirna = lines[i + 2].split(":", 1)[1].strip()
                mfe = lines[i + 5].split(":", 1)[1].strip().split()[0]
                pval = float(lines[i + 6].split(":", 1)[1].strip())
                start = int(lines[i + 8].split()[1]) + 1
                aln = lines[i + 10].strip()
                end = start + len(aln) - 1
                if pval < pvalue_cutoff:
                    fout.write(f"{mirna}\t{target}\t{start}\t{end}\t{mfe}\t{pval}\n")
            except (IndexError, ValueError):
                continue


def read_fasta(path):
    seqs = {}
    cur = None
    parts = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(parts).upper()
                cur = line[1:].split()[0]
                parts = []
            elif cur:
                parts.append(line)
        if cur:
            seqs[cur] = "".join(parts).upper()
    return seqs


def run_targetscan(seed_file, lncrna_fa, out):
    """Simplified TargetScan: match reverse-complement of each miRNA seed
    (2-8 nt) against each lncRNA sequence."""
    seeds = {}
    with open(seed_file) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                seeds[p[0]] = p[1].upper()

    targets = read_fasta(lncrna_fa)
    with open(out, "w") as fout:
        for target, seq in targets.items():
            for mirna, seed in seeds.items():
                motif = revcomp(seed)
                pos = seq.find(motif)
                while pos != -1:
                    fout.write(f"{target}\t{mirna}\t{pos + 1}\t{pos + len(motif)}\n")
                    pos = seq.find(motif, pos + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mature-fa", required=True)
    ap.add_argument("--lncrna-fa", required=True)
    ap.add_argument("--seed-file", required=True)
    ap.add_argument("--miranda-bin", default="miranda")
    ap.add_argument("--rnahybrid-bin", default="RNAhybrid")
    ap.add_argument("--miranda-score", type=float, default=140.0)
    ap.add_argument("--miranda-energy", type=float, default=-10.0)
    ap.add_argument("--rnahybrid-pvalue", type=float, default=0.05)
    ap.add_argument("--rnahybrid-model", default="3utr_human")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    miranda_dir = os.path.join(args.output_dir, "miranda")
    rnahybrid_dir = os.path.join(args.output_dir, "rnahybrid")
    targetscan_dir = os.path.join(args.output_dir, "targetscan")
    for d in (miranda_dir, rnahybrid_dir, targetscan_dir):
        os.makedirs(d, exist_ok=True)

    miranda_raw = os.path.join(miranda_dir, "miranda_interaction_out.txt")
    run([args.miranda_bin, args.mature_fa, args.lncrna_fa,
         "-sc", str(args.miranda_score), "-en", str(args.miranda_energy),
         "-out", miranda_raw])
    parse_miranda(miranda_raw,
                  os.path.join(miranda_dir, "miRanda_interaction_out_select_v2.txt"))

    rnahybrid_raw = os.path.join(rnahybrid_dir, "rnahybrid_raw.txt")
    run([args.rnahybrid_bin, "-s", args.rnahybrid_model,
         args.mature_fa, args.lncrna_fa], out_path=rnahybrid_raw)
    parse_rnahybrid(rnahybrid_raw,
                    os.path.join(rnahybrid_dir, "RNAhybrid_out.txt"),
                    args.rnahybrid_pvalue)

    run_targetscan(args.seed_file, args.lncrna_fa,
                   os.path.join(targetscan_dir, "TargetScan_combine_v2.txt"))

    print(f"done: predictions -> {args.output_dir}")


if __name__ == "__main__":
    main()
