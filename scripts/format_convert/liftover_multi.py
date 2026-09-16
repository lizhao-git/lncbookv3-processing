#!/usr/bin/env python3
"""Lift many BED/bedGraph datasets across assemblies/species from one manifest.

Manifest schema (JSON):
{
  "min_match": 0.95,
  "chains": [
    {"species": "human", "from": "hg19", "to": "hg38", "chain": "chains/hg19ToHg38.over.chain.gz"},
    {"species": "mouse", "from": "mm39", "to": "mm10", "chain": "chains/mm39ToMm10.over.chain.gz"}
  ],
  "datasets": [
    {"name": "human_lncrna", "species": "human", "from": "hg19", "to": "hg38", "input": "human/lncrna.bed"},
    {"name": "mouse_cage",   "species": "mouse", "from": "mm39", "to": "mm10", "input": "mouse/cage.bedgraph",
     "chain": "chains/custom_mm39ToMm10.chain.gz", "min_match": 0.9}
  ]
}
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from genomic_intervals.liftover import run_liftover  # noqa: E402

REPORT_COLUMNS = [
    "dataset", "species", "from", "to", "chain", "status",
    "input_records", "mapped_records", "unmapped_records", "mapping_rate", "message",
]


def _load_manifest(path):
    with open(path, "r", encoding="utf-8") as fh:
        manifest = json.load(fh)
    if not isinstance(manifest, dict):
        raise SystemExit("Manifest must be a JSON object")
    if not manifest.get("datasets"):
        raise SystemExit("Manifest is missing the 'datasets' list")
    return manifest


def _resolve(path, data_root):
    if not path:
        return path
    if data_root and not os.path.isabs(path):
        return os.path.join(data_root, path)
    return path


def _chain_registry(manifest, data_root):
    registry = {}
    for entry in manifest.get("chains", []):
        key = (entry.get("species"), entry.get("from"), entry.get("to"))
        registry[key] = _resolve(entry.get("chain"), data_root)
    return registry


def _safe_name(name, fallback):
    clean = "".join(char if char.isalnum() or char in "._-" else "_" for char in (name or ""))
    return clean or fallback


def main():
    parser = argparse.ArgumentParser(description="Multi-species liftover driver driven by a JSON manifest")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--data-root", default=None)
    parser.add_argument("--output-dir", default="liftover_outputs")
    parser.add_argument("--report", default="liftover_report.tsv")
    args = parser.parse_args()

    manifest = _load_manifest(args.manifest)
    registry = _chain_registry(manifest, args.data_root)
    default_min_match = manifest.get("min_match", 0.95)
    os.makedirs(args.output_dir, exist_ok=True)

    rows = []
    failures = 0
    for index, dataset in enumerate(manifest["datasets"], start=1):
        name = _safe_name(dataset.get("name"), f"dataset{index}")
        species = dataset.get("species", "NA")
        from_build = dataset.get("from", "NA")
        to_build = dataset.get("to", "NA")
        min_match = dataset.get("min_match", default_min_match)
        input_track = _resolve(dataset.get("input"), args.data_root)
        chain = _resolve(dataset.get("chain"), args.data_root) or registry.get((species, from_build, to_build))

        row = {
            "dataset": name, "species": species, "from": from_build, "to": to_build,
            "chain": chain or "NA", "status": "ok", "input_records": "", "mapped_records": "",
            "unmapped_records": "", "mapping_rate": "", "message": "",
        }

        if not input_track:
            row.update(status="failed", message="Dataset is missing 'input'")
        elif not chain:
            row.update(status="failed", message="No chain registered for this species/from/to")
        elif not os.path.exists(input_track):
            row.update(status="failed", message=f"Input does not exist: {input_track}")
        elif not os.path.exists(chain):
            row.update(status="failed", message=f"Chain does not exist: {chain}")
        else:
            try:
                metrics = run_liftover(
                    input_track, chain,
                    os.path.join(args.output_dir, f"{name}.lifted.bed"),
                    os.path.join(args.output_dir, f"{name}.unmapped.bed"),
                    min_match=min_match,
                )
                row.update(metrics)
            except SystemExit as exc:
                row.update(status="failed", message=str(exc))
            except Exception as exc:  # noqa: BLE001 - keep processing the remaining datasets
                row.update(status="failed", message=f"{type(exc).__name__}: {exc}")

        if row["status"] != "ok":
            failures += 1
        rows.append(row)
        print(f"[{index}/{len(manifest['datasets'])}] {name}: {row['status']}" + (f" - {row['message']}" if row["message"] else ""))

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("\t".join(REPORT_COLUMNS) + "\n")
        for row in rows:
            rep.write("\t".join(str(row.get(col, "")) for col in REPORT_COLUMNS) + "\n")

    if failures:
        raise SystemExit(f"{failures} of {len(rows)} dataset(s) failed. See {args.report}")
    print(f"All {len(rows)} dataset(s) lifted successfully. Report: {args.report}")


if __name__ == "__main__":
    main()
