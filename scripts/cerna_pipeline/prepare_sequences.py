#!/usr/bin/env python3
"""Prepare the two sequence inputs for the ceRNA target-prediction tools.

1. ``--gtf`` + ``--genome-fa`` -> ``lncRNA.fa``
   Extract each lncRNA transcript's spliced sequence: exons joined in
   transcript order, reverse-complemented on the ``-`` strand.

2. ``--mature-fa`` -> ``mature_hsa.fa`` + ``mature_seed.txt``
   Keep only human (``Homo sapiens``) mature miRNAs from a miRBase
   ``mature.fa``. U is normalised to T so the query and target alphabets
   match. The seed table (miRNA 2-8 nt, 7-mer) feeds the simplified
   TargetScan step.
"""

import argparse
import os
import re

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def parse_attrs(text):
    attrs = {}
    for m in re.finditer(r'([^\s;]+)\s+"([^"]*)"\s*;', text):
        attrs[m.group(1)] = m.group(2)
    return attrs


def load_fai(fa_path):
    """Return ``{chrom: (base_offset, bases_per_line, line_bytes, length)}``.

    Uses ``<fa>.fai`` when present (samtools format), otherwise scans the
    file once recording offsets only (no sequence is kept in memory).
    """
    fai = {}
    fai_path = fa_path + ".fai"
    if os.path.exists(fai_path):
        with open(fai_path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                fai[p[0]] = (int(p[2]), int(p[3]), int(p[4]), int(p[1]))
        return fai

    fai = {}
    byte = 0
    cur = None
    cur_offset = None
    cur_bpl = None
    cur_lb = None
    cur_len = 0
    with open(fa_path) as fh:
        for raw in fh:
            if raw.startswith(">"):
                if cur is not None:
                    fai[cur] = (cur_offset, cur_bpl, cur_lb, cur_len)
                cur = raw[1:].split()[0]
                cur_offset = byte + len(raw)
                cur_bpl = None
                cur_lb = None
                cur_len = 0
            else:
                line = raw.rstrip("\n")
                if cur_bpl is None:
                    cur_bpl = len(line)
                    cur_lb = len(raw)  # includes the newline
                cur_len += len(line)
            byte += len(raw)
        if cur is not None:
            fai[cur] = (cur_offset, cur_bpl, cur_lb, cur_len)
    return fai


def _match_chrom(chrom, fai):
    if chrom in fai:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in fai:
        return chrom[3:]
    if ("chr" + chrom) in fai:
        return "chr" + chrom
    return None


def fetch_seq(fh, fai, chrom, start_0, end_0):
    """Return uppercase bases for the 0-based half-open interval [start_0, end_0)."""
    offset, bpl, lb, length = fai[chrom]
    s = start_0 // bpl
    e = (end_0 - 1) // bpl
    rel_start = offset + s * lb + (start_0 % bpl)
    n_bytes = (e - s) * lb + ((end_0 - 1) % bpl) + 1 - (start_0 % bpl)
    fh.seek(rel_start)
    return fh.read(n_bytes).replace("\n", "").upper()


def gtf_to_transcripts(gtf, genome_fa, out_fa):
    tx_exons = {}
    with open(gtf) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9 or cols[2] != "exon":
                continue
            attrs = parse_attrs(cols[8])
            tx = attrs.get("transcript_id")
            if not tx:
                continue
            try:
                s = int(cols[3])
                e = int(cols[4])
            except ValueError:
                continue
            tx_exons.setdefault(tx, []).append((cols[0], s, e, cols[6]))

    fai = load_fai(genome_fa)
    with open(genome_fa) as gfh, open(out_fa, "w") as out:
        for tx in sorted(tx_exons):
            exons = sorted(set(tx_exons[tx]), key=lambda x: (x[1], x[2]))
            strand = exons[0][3]
            if strand == "-":
                exons.sort(key=lambda x: (-x[1], -x[2]))
            parts = []
            for (chrom, s, e, _st) in exons:
                chrom = _match_chrom(chrom, fai)
                if chrom is None:
                    continue
                sub = fetch_seq(gfh, fai, chrom, s - 1, e)
                parts.append(revcomp(sub) if strand == "-" else sub)
            out.write(f">{tx}\n{''.join(parts)}\n")


def mature_to_hsa(mature_fa, out_fa, out_seed):
    seen = set()
    cur_id = None
    seqs = []

    def flush():
        nonlocal cur_id, seqs
        if cur_id and cur_id not in seen:
            seen.add(cur_id)
            seq = "".join(seqs).upper().replace("U", "T")
            if len(seq) >= 8:
                out.write(f">{cur_id}\n{seq}\n")
                seed_fh.write(f"{cur_id}\t{seq[1:8]}\t9606\n")
        cur_id = None
        seqs = []

    with open(mature_fa) as fin, open(out_fa, "w") as out, open(out_seed, "w") as seed_fh:
        for raw in fin:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                if "Homo sapiens" in header:
                    cur_id = header.split()[0]
            elif cur_id:
                seqs.append(line)
        flush()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf")
    ap.add_argument("--genome-fa")
    ap.add_argument("--mature-fa")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    if args.gtf and args.genome_fa:
        gtf_to_transcripts(args.gtf, args.genome_fa,
                           os.path.join(args.output_dir, "lncRNA.fa"))
    if args.mature_fa:
        mature_to_hsa(args.mature_fa,
                      os.path.join(args.output_dir, "mature_hsa.fa"),
                      os.path.join(args.output_dir, "mature_seed.txt"))
    print(f"done: prepared sequences -> {args.output_dir}")


if __name__ == "__main__":
    main()
