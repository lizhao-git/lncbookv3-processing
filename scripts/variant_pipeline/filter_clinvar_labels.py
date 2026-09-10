#!/usr/bin/env python3
import argparse
from collections import Counter


TARGET_LABELS = {
    "benign": "Benign",
    "pathogenic": "Pathogenic",
    "drug response": "Drug response",
    "protective": "Protective",
    "risk factor": "Risk factor",
    "affects": "Affects",
    "affectes": "Affects",
}

REFSEQ_ACCESSION = {
    "chr1": "NC_000001.11",
    "chr2": "NC_000002.12",
    "chr3": "NC_000003.12",
    "chr4": "NC_000004.12",
    "chr5": "NC_000005.10",
    "chr6": "NC_000006.12",
    "chr7": "NC_000007.14",
    "chr8": "NC_000008.11",
    "chr9": "NC_000009.12",
    "chr10": "NC_000010.11",
    "chr11": "NC_000011.10",
    "chr12": "NC_000012.12",
    "chr13": "NC_000013.11",
    "chr14": "NC_000014.9",
    "chr15": "NC_000015.10",
    "chr16": "NC_000016.10",
    "chr17": "NC_000017.11",
    "chr18": "NC_000018.10",
    "chr19": "NC_000019.10",
    "chr20": "NC_000020.11",
    "chr21": "NC_000021.9",
    "chr22": "NC_000022.11",
    "chrX": "NC_000023.11",
    "chrY": "NC_000024.10",
    "chrM": "NC_012920.1",
}


def normalize_chrom(chrom: str) -> str:
    c = chrom.strip()
    if c.startswith("chr"):
        return c
    if c == "MT":
        return "chrM"
    return f"chr{c}"


def parse_info(info: str):
    items = {}
    for token in info.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            items[k] = v
        else:
            items[token] = True
    return items


def pick_labels(clnsig: str):
    s = clnsig.lower().replace("_", " ")
    found = []
    for k, canonical in TARGET_LABELS.items():
        if k in s and canonical not in found:
            found.append(canonical)
    return found


def build_variant_name(chrom_norm: str, pos: int, ref: str, alt: str) -> str:
    refseq = REFSEQ_ACCESSION.get(chrom_norm)
    alt_first = alt.split(",", 1)[0].strip() if alt else ""
    alt_norm = alt_first if alt_first and alt_first != "." else "?"
    if refseq:
        return f"{refseq}:g.{pos}{ref}>{alt_norm}"
    return f"{chrom_norm}:g.{pos}{ref}>{alt_norm}"


def main():
    parser = argparse.ArgumentParser(description="Filter ClinVar VCF by definite clinical significance labels")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    kept = 0
    counts = Counter()

    with open(args.input_vcf, "r", encoding="utf-8") as inp, \
         open(args.output_bed, "w", encoding="utf-8") as bed, \
         open(args.output_vcf, "w", encoding="utf-8") as vcf:

        for raw in inp:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                vcf.write(raw)
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
            chrom_norm = normalize_chrom(chrom)
            info_map = parse_info(info)
            clnsig = str(info_map.get("CLNSIG", ""))
            labels = pick_labels(clnsig)
            if not labels:
                continue

            allele_id = str(info_map.get("ALLELEID", "NA"))
            disease_name = str(info_map.get("CLNDN", "NA"))
            dbsnp_from_info = str(info_map.get("RS", "")).strip()

            try:
                pos_i = int(pos)
            except ValueError:
                continue

            start0 = pos_i - 1
            end = start0 + max(1, len(ref))
            variant_id = vid if vid and vid != "." else f"{chrom_norm}:{pos}:{ref}:{alt}"
            dbsnp_id = vid if vid.lower().startswith("rs") else (f"rs{dbsnp_from_info}" if dbsnp_from_info else "NA")
            variant_name = build_variant_name(chrom_norm, pos_i, ref, alt)
            label_text = "|".join(labels)

            bed.write(
                f"{chrom_norm}\t{start0}\t{end}\t{variant_id}\t{label_text}\t{allele_id}\t{disease_name}\t{variant_name}\t{dbsnp_id}\n"
            )
            vcf.write(raw)

            kept += 1
            for label in labels:
                counts[label] += 1

    with open(args.report, "w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"kept_variants\t{kept}\n")
        for label in ["Pathogenic", "Benign", "Affects", "Drug response", "Protective", "Risk factor"]:
            rep.write(f"label_{label.replace(' ', '_')}\t{counts.get(label, 0)}\n")


if __name__ == "__main__":
    main()
