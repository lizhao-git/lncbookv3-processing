#!/usr/bin/env python3
"""ceRNA pipeline orchestrator (Nextflow/Docker entry point).

Reads a single JSON manifest and drives the refactored ceRNA stages:

    parse_tools -> predict_ceRNA -> experiment_annotate

Manifest schema::

    {
      "miranda":               "path/miRanda_interaction_out_select_v2.txt",
      "targetscan":            "path/TargetScan_combine_v2.txt",
      "rnahybrid":             "path/RNAhybrid_split",   // dir of *_new_out.txt or file
      "transcript_id_alias":   "path/transcript_id_alias.txt",
      "hgnc_symbol_transcript": "path/HGNC_symbol_transcript.txt",
      "transcript_gene":       "path/transcript_gene.txt",
      "lncrnawiki":            "path/LncRNAWiki2.0-ceRNA.txt",
      "hmdd":                  "path/alldata_v2.txt"
    }

All manifest paths are resolved against ``--data-root`` (a workflow data-root mount)
when relative; absolute paths are used as-is.
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
    ap.add_argument("--output-dir", default="cerna_outputs")
    args = ap.parse_args()

    with open(args.manifest) as fh:
        manifest = json.load(fh)

    data_root = args.data_root
    manifest_dir = os.path.dirname(os.path.abspath(args.manifest))
    out_root = args.output_dir
    os.makedirs(out_root, exist_ok=True)

    miranda = _resolve(manifest.get("miranda"), data_root, manifest_dir)
    targetscan = _resolve(manifest.get("targetscan"), data_root, manifest_dir)
    rnahybrid = _resolve(manifest.get("rnahybrid"), data_root, manifest_dir)

    # 1) normalize the three tools' outputs
    norm_dir = os.path.join(out_root, "normalized")
    _python(os.path.join(MODULE_DIR, "parse_tools.py"),
            "--miranda", miranda,
            "--targetscan", targetscan,
            "--rnahybrid", rnahybrid,
            "--output-dir", norm_dir)

    # 2) 3-way intersection + binding-site refinement
    predict_dir = os.path.join(out_root, "predict")
    _python(os.path.join(MODULE_DIR, "predict_ceRNA.py"),
            "--miranda", os.path.join(norm_dir, "miranda.tsv"),
            "--targetscan", os.path.join(norm_dir, "targetscan.tsv"),
            "--rnahybrid", os.path.join(norm_dir, "rnahybrid.tsv"),
            "--output-dir", predict_dir)

    # 3) experimental validation + final table + disease annotation
    annotate_dir = os.path.join(out_root, "annotate")
    _python(os.path.join(MODULE_DIR, "experiment_annotate.py"),
            "--lncrnawiki", _resolve(manifest["lncrnawiki"], data_root, manifest_dir),
            "--hgnc", _resolve(manifest["hgnc_symbol_transcript"], data_root, manifest_dir),
            "--transcript-id-alias", _resolve(manifest["transcript_id_alias"], data_root, manifest_dir),
            "--transcript-gene", _resolve(manifest["transcript_gene"], data_root, manifest_dir),
            "--hmdd", _resolve(manifest["hmdd"], data_root, manifest_dir),
            "--three-overlap", os.path.join(predict_dir, "three_overlap_unique.txt"),
            "--output-dir", annotate_dir)

    print(f"done: ceRNA outputs -> {out_root}", file=sys.stderr)


if __name__ == "__main__":
    main()
