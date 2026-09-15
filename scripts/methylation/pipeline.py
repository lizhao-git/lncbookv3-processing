#!/usr/bin/env python3
"""Methylation pipeline orchestrator (Nextflow/Docker entry point).

Reads a single JSON manifest describing every GEO/TCGA dataset plus the
cross-dataset aggregation config, then drives the refactored stages in order:

    per dataset:  preprocess -> extract -> [beta_convert] -> difftest (body + promoter)
    then:         aggregate  (cross-dataset significance matrix + labels)

This replaces the legacy ``run_legacy_analyze_and_format.py`` wrapper: it runs
the refactored, path-parameterized Python 3 modules instead of staging a
hard-coded ``/disk1/...`` tree and executing Python 2 scripts.

Manifest schema::

    {
      "gene_info":     "path/to/gene_info.txt",           // 1-based closed
      "promoter_info": "path/to/gene_promoter_info.txt",  // optional
      "gene_list":     "path/to/general_gene_list.txt",   // optional (GTF-derived)
      "datasets": [
        {
          "name": "GSE109875",
          "disease": "Autism",
          "category": "brain",               // brain | cancer
          "format": "bismark",               // bigwig | bigbed | bismark | bed
          "input_dir": "path/to/input",      // per-dataset raw input
          "chain": null,                     // optional hg19->hg38 chain file
          "value_index": null,               // optional value column override
          "class_map": {"GSMx": "case", ...},  // optional; OR class_regex below
          "class_regex": {"case": "MB|tumor", "control": "Normal"},
          "method": "wilcox",                // wilcox | foldchange
          "case_label": "case",
          "control_label": "control",
          "pval": 0.02,
          "min_median": 0.2,
          "fc": 2.0,
          "min_abs": 20.0,
          "beta_scale": null                 // optional, e.g. 0.01 (percent->beta)
        }
      ]
    }

All manifest paths are resolved against ``--data-root`` (a workflow data-root mount)
when they are relative; absolute paths are used as-is.
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


def _class_args(ds, workdir):
    if ds.get("class_map"):
        cmap_file = os.path.join(workdir, f"{ds['name']}_class_map.json")
        with open(cmap_file, "w") as fh:
            json.dump(ds["class_map"], fh)
        return ["--class-map", cmap_file]
    if ds.get("class_regex"):
        specs = ";".join(f"{cls}={pattern}"
                         for cls, pattern in ds["class_regex"].items())
        return ["--class-regex", specs]
    return []


def run_dataset(ds, gene_info, promoter_info, out_root, data_root, manifest_dir):
    name = ds["name"]
    workdir = os.path.join(out_root, name)
    os.makedirs(workdir, exist_ok=True)

    input_dir = _resolve(ds["input_dir"], data_root, manifest_dir)
    chain = _resolve(ds.get("chain"), data_root, manifest_dir)

    # 1) preprocess
    pre_dir = os.path.join(workdir, "pre")
    pre_args = ["--format", ds["format"], "--input-dir", input_dir,
                "--output-dir", pre_dir]
    if ds.get("value_index") is not None:
        pre_args += ["--value-index", str(ds["value_index"])]
    if chain:
        pre_args += ["--chain", chain]
    _python(os.path.join(MODULE_DIR, "preprocess.py"), *pre_args)

    bed_dir = os.path.join(pre_dir, "hg38_output" if chain else "bed")

    # 2) extract (gene body + promoter matrices)
    ext_args = ["--bed-dir", bed_dir, "--gene-info", gene_info,
                "--prefix", name, "--output-dir", workdir]
    if promoter_info:
        ext_args += ["--promoter-info", promoter_info]
    ext_args += _class_args(ds, workdir)
    _python(os.path.join(MODULE_DIR, "extract.py"), *ext_args)

    body_matrix = os.path.join(workdir, f"{name}_body")
    promoter_matrix = os.path.join(workdir, f"{name}_promoter")

    # 3) optional beta rescale (before differential testing)
    if ds.get("beta_scale"):
        scale = ds["beta_scale"]
        _python(os.path.join(MODULE_DIR, "beta_convert.py"),
                "--input", body_matrix, "--output", body_matrix + "_0_1", "--scale", str(scale))
        _python(os.path.join(MODULE_DIR, "beta_convert.py"),
                "--input", promoter_matrix, "--output", promoter_matrix + "_0_1", "--scale", str(scale))
        body_matrix += "_0_1"
        promoter_matrix += "_0_1"

    # 4) differential testing (body + promoter)
    common_dt = ["--prefix", name, "--output-dir", workdir,
                 "--method", ds.get("method", "wilcox"),
                 "--case-label", ds.get("case_label", "case"),
                 "--control-label", ds.get("control_label", "control")]
    if ds.get("method", "wilcox") == "wilcox":
        common_dt += ["--pval", str(ds.get("pval", 0.02)),
                      "--min-median", str(ds.get("min_median", 0.2))]
    else:
        common_dt += ["--fc", str(ds.get("fc", 2.0)),
                      "--min-abs", str(ds.get("min_abs", 20.0))]

    _python(os.path.join(MODULE_DIR, "difftest.py"),
            "--matrix", body_matrix, "--region", "body", *common_dt)
    _python(os.path.join(MODULE_DIR, "difftest.py"),
            "--matrix", promoter_matrix, "--region", "promoter", *common_dt)

    return {
        "name": name,
        "disease": ds["disease"],
        "dataset": ds.get("dataset", name),
        "category": ds.get("category", "cancer"),
        "body_sig": os.path.join(workdir, f"{name}_body_sig.txt"),
        "promoter_sig": os.path.join(workdir, f"{name}_promoter_sig.txt"),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--data-root",
                    help="base directory for resolving relative manifest paths")
    ap.add_argument("--output-dir", default="methylation_outputs")
    ap.add_argument("--report", default="methylation_report.tsv")
    args = ap.parse_args()

    with open(args.manifest) as fh:
        manifest = json.load(fh)

    data_root = args.data_root
    manifest_dir = os.path.dirname(os.path.abspath(args.manifest))
    out_root = args.output_dir
    os.makedirs(out_root, exist_ok=True)

    gene_info = _resolve(manifest["gene_info"], data_root, manifest_dir)
    promoter_info = _resolve(manifest.get("promoter_info"), data_root, manifest_dir)
    gene_list = _resolve(manifest.get("gene_list"), data_root, manifest_dir)

    report_rows = []
    results = []
    for ds in manifest["datasets"]:
        status = "ok"
        try:
            res = run_dataset(ds, gene_info, promoter_info, out_root, data_root, manifest_dir)
            results.append(res)
        except subprocess.CalledProcessError as exc:
            status = f"failed: {exc}"
        report_rows.append((ds["name"], ds["disease"], status))

    # 5) cross-dataset aggregation
    agg_manifest = {
        "gene_list": gene_list,
        "diseases": [
            {
                "name": r["disease"],
                "dataset": r["dataset"],
                "category": r["category"],
                "body": r["body_sig"],
                "promoter": r["promoter_sig"],
            }
            for r in results
        ],
    }
    agg_manifest_path = os.path.join(out_root, "aggregate_manifest.json")
    with open(agg_manifest_path, "w") as fh:
        json.dump(agg_manifest, fh, indent=2)

    agg_dir = os.path.join(out_root, "aggregate")
    _python(os.path.join(MODULE_DIR, "aggregate.py"),
            "--manifest", agg_manifest_path, "--output-dir", agg_dir)

    with open(args.report, "w") as fh:
        fh.write("dataset\tdisease\tstatus\n")
        for row in report_rows:
            fh.write("\t".join(row) + "\n")

    failed = [row for row in report_rows if row[2] != "ok"]
    print(f"done: {len(results)} datasets, aggregate -> {agg_dir}", file=sys.stderr)
    if failed:
        raise SystemExit(f"{len(failed)} dataset(s) failed: {failed[0][0]}")


if __name__ == "__main__":
    main()
