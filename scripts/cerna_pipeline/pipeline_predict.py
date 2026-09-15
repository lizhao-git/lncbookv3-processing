#!/usr/bin/env python3
"""ceRNA target-prediction orchestrator (Nextflow/Docker entry point).

Reads a single JSON manifest and drives the upstream prediction stage:

    prepare_sequences -> predict_targets

The resulting three normalized files (under <output-dir>/predictions/) are
the inputs referenced by the downstream ceRNA pipeline manifest
(``cerna.example.json``).

Manifest schema::

    {
      "genome_fa":        "path/hg38.fa",
      "gtf":              "path/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf",
      "mature_mirna_fa":  "path/mature.fa",
      "miranda_bin":      "miranda",
      "rnahybrid_bin":    "RNAhybrid",
      "miranda_score":     140.0,
      "miranda_energy":    -10.0,
      "rnahybrid_pvalue":  0.05,
      "rnahybrid_model":   "3utr_human"
    }

Relative manifest paths are resolved against ``--data-root`` (a workflow
Directory mount); absolute paths are used as-is.
"""

import argparse
import json
import os
import subprocess
import sys

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


def _python(*argv):
    cmd = [sys.executable] + list(argv)
    print("+ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def _resolve(path, data_root, manifest_dir):
    if not path:
        return None
    if os.path.isabs(path):
        return path
    for base in (data_root, manifest_dir):
        if base:
            candidate = os.path.join(base, path)
            if os.path.exists(candidate):
                return candidate
    return os.path.join(data_root or manifest_dir or ".", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--data-root",
                    help="base directory for resolving relative manifest paths")
    ap.add_argument("--output-dir", default="cerna_predict")
    args = ap.parse_args()

    with open(args.manifest) as fh:
        manifest = json.load(fh)

    data_root = args.data_root
    manifest_dir = os.path.dirname(os.path.abspath(args.manifest))
    out_root = args.output_dir
    os.makedirs(out_root, exist_ok=True)

    genome = _resolve(manifest.get("genome_fa"), data_root, manifest_dir)
    gtf = _resolve(manifest.get("gtf"), data_root, manifest_dir)
    mature = _resolve(manifest.get("mature_mirna_fa"), data_root, manifest_dir)

    # 1) prepare sequences
    seq_dir = os.path.join(out_root, "sequences")
    prep = [os.path.join(MODULE_DIR, "prepare_sequences.py"),
            "--output-dir", seq_dir]
    if gtf and genome:
        prep += ["--gtf", gtf, "--genome-fa", genome]
    if mature:
        prep += ["--mature-fa", mature]
    _python(*prep)

    # 2) run the three tools
    pred_dir = os.path.join(out_root, "predictions")
    _python(os.path.join(MODULE_DIR, "predict_targets.py"),
            "--mature-fa", os.path.join(seq_dir, "mature_hsa.fa"),
            "--lncrna-fa", os.path.join(seq_dir, "lncRNA.fa"),
            "--seed-file", os.path.join(seq_dir, "mature_seed.txt"),
            "--miranda-bin", manifest.get("miranda_bin", "miranda"),
            "--rnahybrid-bin", manifest.get("rnahybrid_bin", "RNAhybrid"),
            "--miranda-score", str(manifest.get("miranda_score", 140.0)),
            "--miranda-energy", str(manifest.get("miranda_energy", -10.0)),
            "--rnahybrid-pvalue", str(manifest.get("rnahybrid_pvalue", 0.05)),
            "--rnahybrid-model", manifest.get("rnahybrid_model", "3utr_human"),
            "--output-dir", pred_dir)

    print(f"done: ceRNA predictions -> {out_root}", file=sys.stderr)


if __name__ == "__main__":
    main()
