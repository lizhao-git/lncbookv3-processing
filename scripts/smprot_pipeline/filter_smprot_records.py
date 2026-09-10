#!/usr/bin/env python3
import argparse
import csv

COLUMN_ALIASES = {
    "chrom": ["chrom", "chr", "chromosome"],
    "start": ["start", "chromstart", "tx_start", "pos", "begin"],
    "end": ["end", "chromend", "tx_end", "stop"],
    "protein_id": ["protein_id", "smprot_id", "smprotid", "proteinid", "id", "protein"],
    "evidence": ["evidence", "evidence_type", "support", "support_evidence"],
    "source": ["source", "dataset", "source_db"],
}


def normalize(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_")


def resolve_columns(fieldnames):
    normalized = {normalize(col): col for col in fieldnames}
    resolved = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            key = normalize(alias)
            if key in normalized:
                resolved[canonical] = normalized[key]
                break
    return resolved


def normalize_chrom(chrom: str) -> str:
    c = chrom.strip()
    if c.startswith("chr"):
        return c
    if c == "MT":
        return "chrM"
    return f"chr{c}"


def main():
    parser = argparse.ArgumentParser(description="Filter/normalize SmProt coordinates and emit BED")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    kept = 0
    skipped = 0

    with open(args.input_tsv, "r", encoding="utf-8") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        resolved = resolve_columns(reader.fieldnames or [])
        required = ["chrom", "start", "end", "protein_id"]
        missing = [name for name in required if name not in resolved]
        if missing:
            raise SystemExit(f"Missing required SmProt columns: {', '.join(missing)}")

        with open(args.output_bed, "w", encoding="utf-8") as bed, open(args.output_tsv, "w", encoding="utf-8", newline="") as out_tsv:
            writer = csv.DictWriter(out_tsv, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()

            for row in reader:
                chrom_raw = str(row[resolved["chrom"]]).strip()
                start_raw = str(row[resolved["start"]]).strip()
                end_raw = str(row[resolved["end"]]).strip()
                protein_id = str(row[resolved["protein_id"]]).strip()
                if not chrom_raw or not start_raw or not end_raw:
                    skipped += 1
                    continue
                try:
                    start = int(start_raw)
                    end = int(end_raw)
                except ValueError:
                    skipped += 1
                    continue
                if start < 0 or end <= start:
                    skipped += 1
                    continue
                if not protein_id:
                    skipped += 1
                    continue

                evidence = str(row[resolved["evidence"]]).strip() if "evidence" in resolved else "NA"
                source = str(row[resolved["source"]]).strip() if "source" in resolved else "SmProt"
                chrom = normalize_chrom(chrom_raw)

                bed.write(f"{chrom}\t{start}\t{end}\t{protein_id}\t{evidence or 'NA'}\t{source or 'SmProt'}\n")
                writer.writerow(row)
                kept += 1

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"kept_records\t{kept}\n")
        rep.write(f"skipped_records\t{skipped}\n")


if __name__ == "__main__":
    main()
