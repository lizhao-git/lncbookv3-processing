#!/usr/bin/env python3
"""Conservation pipeline orchestrator (Nextflow/Docker entry point).

Manifest schema::

    {
      "thresholds": {
        "min_mapped_length": 50,
        "min_transcript_coverage": 0.2,
        "min_intron_coverage": 0.0
      },
      "species_order": "species_order.tsv",
      "synteny_anchors": "synteny_anchors.tsv",
      "species": [
        {
          "prefix": "mouse",
          "species": "Mus_musculus",
          "lnc": {
            "query_fasta": "fas/mouse_lnc_query_simply.fa",
            "target_fasta": "fas/mouse_lnc_target_simply.fa",
            "group_info": "lnc_group_info/mouse_lnc_group_info.txt",
            "lengths": "lnc_lengths/mouse_lnc_lengths.txt",
            "exon_positions": "exon_positions.txt",
            "target_info": "bedtool_exonCDS/mouse_lnc_all_target_info.txt",
            "target_addition": "CDSexon_range/mouse_exoncount_ORFexonrange.txt",
            "transcript_lengths": "all_lnc_length.txt",
            "transcript_meta": "hg38_all_trans_meta_info_new.txt"
          },
          "pc": {
            "...": "optional protein-coding control inputs"
          }
        }
      ]
    }
"""

import argparse
import os
import sys
from argparse import Namespace

if __package__ is None or __package__ == "":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from conservation.classify import classify_rows
from conservation.common import ensure_dir, load_manifest, read_tsv, resolve_path, write_tsv
from conservation.gene_age import assign_gene_age
from conservation.homology import build_homology
from conservation.postprocess import run_postprocess
from conservation.synteny import call_synteny


def resolved_config(config, data_root, manifest_dir):
    if not config:
        return {}
    out = {}
    path_keys = {
        "query_fasta", "target_fasta", "group_info", "lengths", "exon_positions",
        "target_info", "target_addition", "transcript_lengths", "transcript_meta",
        "cds_positions",
    }
    for key, value in config.items():
        out[key] = resolve_path(value, data_root, manifest_dir) if key in path_keys else value
    return out


def postprocess_one(mode, prefix, species, config, out_root):
    species_dir = ensure_dir(os.path.join(out_root, "statistics", species))
    output = os.path.join(species_dir, f"{prefix}_{mode}_all_statistics.tsv")
    report = os.path.join(species_dir, f"{prefix}_{mode}_postprocess_report.tsv")
    args = Namespace(
        mode=mode,
        query_fasta=config["query_fasta"],
        target_fasta=config["target_fasta"],
        group_info=config["group_info"],
        lengths=config["lengths"],
        exon_positions=config["exon_positions"],
        target_info=config.get("target_info"),
        target_addition=config.get("target_addition"),
        transcript_lengths=config.get("transcript_lengths"),
        transcript_meta=config.get("transcript_meta"),
        cds_positions=config.get("cds_positions"),
        overlap_cutoff=float(config.get("overlap_cutoff", 0.5)),
        nearby_distance=int(config.get("nearby_distance", 100000)),
        output=output,
        report=report,
    )
    run_postprocess(args)
    return output, report


def combine(files):
    rows = []
    for path in files:
        rows.extend(read_tsv(path))
    return rows


def summarize_qc(report_files, output):
    rows = []
    for path in report_files:
        for row in read_tsv(path):
            row = dict(row)
            row["report_file"] = path
            rows.append(row)
    write_tsv(output, rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--data-root", help="base directory for resolving relative manifest paths")
    ap.add_argument("--output-dir", default="conservation_outputs")
    ap.add_argument("--report", default="conservation_report.tsv")
    args = ap.parse_args()

    manifest = load_manifest(args.manifest)
    manifest_dir = os.path.dirname(os.path.abspath(args.manifest))
    out_root = ensure_dir(args.output_dir)
    thresholds = manifest.get("thresholds", {})
    min_mapped = float(thresholds.get("min_mapped_length", 50))
    min_tx_cov = float(thresholds.get("min_transcript_coverage", 0.2))
    min_intron_cov = float(thresholds.get("min_intron_coverage", 0.0))

    lnc_stats = []
    pc_stats = []
    reports = []
    run_rows = []

    for item in manifest.get("species", []):
        species = item.get("species") or item["prefix"]
        prefix = item.get("prefix") or species
        status = "ok"
        try:
            if item.get("lnc"):
                config = resolved_config(item["lnc"], args.data_root, manifest_dir)
                output, report = postprocess_one("lnc", prefix, species, config, out_root)
                lnc_stats.append(output)
                reports.append(report)
            if item.get("pc"):
                config = resolved_config(item["pc"], args.data_root, manifest_dir)
                output, report = postprocess_one("pc", prefix, species, config, out_root)
                pc_stats.append(output)
                reports.append(report)
        except Exception as exc:
            status = f"failed: {exc}"
        run_rows.append({"prefix": prefix, "species": species, "status": status})

    lnc_rows = []
    for path in lnc_stats:
        species = os.path.basename(os.path.dirname(path))
        lnc_rows.extend(classify_rows(read_tsv(path), species, min_mapped, min_tx_cov, min_intron_cov))

    sequence_out = os.path.join(out_root, "sequence_conservation.tsv")
    preferred = [
        "target_species", "ids", "grouped_ids", "geneID", "gene_name",
        "transcript_length", "true_mapped_length", "mapped_length",
        "matched_length", "match_ratio", "transcript_coverage",
        "exon_covered_ratio", "intron_coverage", "pass_length_filter",
        "pass_transcript_coverage_filter", "pass_intron_coverage_filter",
        "is_conserved",
    ]
    write_tsv(sequence_out, lnc_rows, preferred + sorted({k for r in lnc_rows for k in r if k not in preferred}))

    evidence, best = build_homology(lnc_rows)
    homology_all = os.path.join(out_root, "gene_homology.tsv")
    homology_best = os.path.join(out_root, "gene_homology_best_hit.tsv")
    write_tsv(homology_all, evidence)
    write_tsv(homology_best, best)

    synteny_out = os.path.join(out_root, "synteny.tsv")
    anchors = resolve_path(manifest.get("synteny_anchors"), args.data_root, manifest_dir)
    if anchors and os.path.exists(anchors):
        write_tsv(synteny_out, call_synteny(best, read_tsv(anchors)))
    else:
        write_tsv(synteny_out, [dict(row, syntenic_support="") for row in best])

    species_order = resolve_path(manifest.get("species_order"), args.data_root, manifest_dir)
    gene_age_out = os.path.join(out_root, "gene_age.tsv")
    if species_order and os.path.exists(species_order):
        write_tsv(gene_age_out, assign_gene_age(best, species_order))
    else:
        write_tsv(gene_age_out, [])

    if pc_stats:
        write_tsv(os.path.join(out_root, "pc_control_statistics.tsv"), combine(pc_stats))

    summarize_qc(reports, os.path.join(out_root, "conservation_qc.tsv"))
    write_tsv(args.report, run_rows, ["prefix", "species", "status"])

    failed = [row for row in run_rows if row["status"] != "ok"]
    print(f"done: conservation outputs -> {out_root}", file=sys.stderr)
    if failed:
        raise SystemExit(f"{len(failed)} species failed: {failed[0]['species']}")


if __name__ == "__main__":
    main()

