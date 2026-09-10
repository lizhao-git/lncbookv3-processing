#!/usr/bin/env python3
"""Stage C: cross-dataset aggregation, counting, labelling and the
database-import pivot.

Consumes per-dataset significant-call files produced by ``difftest.py``
(``gene  direction  evidence``) and merges them into a gene x disease matrix.

Manifest (JSON) format::

  {
    "gene_list": "path/to/general_gene_list.txt", // optional, one id per line
    "diseases": [
      {
        "name": "Autism",
        "dataset": "GSE109875",
        "category": "brain",                  // brain | cancer
        "body": "path/GSE109875_body_sig.txt",
        "promoter": "path/GSE109875_promoter_sig.txt"
      }
    ]
  }

Outputs (written to --output-dir):
  * Methylation_sig_body_out.txt / _promoter_out.txt      gene x disease matrix
  * body_disease_type.txt / promoter_disease_type.txt     gene -> category
  * Methylation_sig_body_out_number.txt (+_pandisease)    >=5 & 0 filtering
  * Sig_label_out.txt                                     long-form labels
  * Sig_label_out_database_import.txt                     wide import matrix
"""

import argparse
import json
import os
import sys


def load_sig(path):
    """Return ``{gene: (direction, evidence)}`` from a sig file."""
    result = {}
    if not path or not os.path.exists(path):
        return result
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene = parts[0]
            direction = "High" if parts[1] == "High" else "Low"
            evidence = parts[2] if len(parts) > 2 else ""
            result[gene] = (direction, evidence)
    return result


def load_gene_list(path):
    if not path:
        return None
    ids = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                ids.append(parts[0])
    return set(ids)


def direction_label(direction):
    return "hyper" if direction == "High" else "hypo"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    with open(args.manifest) as fh:
        manifest = json.load(fh)

    gene_list = load_gene_list(manifest.get("gene_list"))
    diseases = manifest["diseases"]
    os.makedirs(args.output_dir, exist_ok=True)

    # body/promoter gene -> {disease: direction}; evidence for labels
    body = {}
    promoter = {}
    for d in diseases:
        name = d["name"]
        body[name] = load_sig(d.get("body"))
        promoter[name] = load_sig(d.get("promoter"))

    def keep(gene):
        if gene_list is not None and gene not in gene_list:
            return False
        return True

    def genes_in(region_dict):
        gs = set()
        for d in diseases:
            gs.update(region_dict[d["name"]].keys())
        return sorted(g for g in gs if keep(g))

    body_genes = genes_in(body)
    promoter_genes = genes_in(promoter)
    all_genes = sorted(set(body_genes) | set(promoter_genes))

    def write_matrix(path, region_dict, genes):
        with open(path, "w") as fh:
            fh.write("Gene\t" + "\t".join(d["name"] for d in diseases) + "\n")
            for g in genes:
                cells = []
                for d in diseases:
                    hit = region_dict[d["name"]].get(g)
                    cells.append(direction_label(hit[0]) if hit else "-")
                fh.write(g + "\t" + "\t".join(cells) + "\n")

    write_matrix(os.path.join(args.output_dir, "Methylation_sig_body_out.txt"),
                 body, body_genes)
    write_matrix(os.path.join(args.output_dir, "Methylation_sig_promoter_out.txt"),
                 promoter, promoter_genes)

    # disease-type annotations (cancer/brain per gene, per region)
    def write_disease_type(path, region_dict, genes, region):
        with open(path, "w") as fh:
            for g in genes:
                cats = set()
                for d in diseases:
                    if g in region_dict[d["name"]]:
                        cats.add(d["category"])
                for cat in sorted(cats):
                    fh.write(f"{g}\t{cat}\t{region}\n")

    write_disease_type(os.path.join(args.output_dir, "body_disease_type.txt"),
                       body, body_genes, "body")
    write_disease_type(os.path.join(args.output_dir, "promoter_disease_type.txt"),
                       promoter, promoter_genes, "promoter")

    # counting / filtering: hyper >=5 & hypo==0 (or hypo >=5 & hyper==0)
    def count_filter(path, region_dict, genes, diseases_subset):
        with open(path, "w") as fh:
            for g in genes:
                h = l = 0
                for d in diseases_subset:
                    hit = region_dict[d["name"]].get(g)
                    if not hit:
                        continue
                    if hit[0] == "High":
                        h += 1
                    else:
                        l += 1
                if h >= 5 and l == 0:
                    fh.write(f"{g}\thyper\n")
                elif l >= 5 and h == 0:
                    fh.write(f"{g}\thypo\n")

    cancer_diseases = [d for d in diseases if d["category"] == "cancer"]
    count_filter(os.path.join(args.output_dir, f"Methylation_sig_body_out_number.txt"),
                 body, body_genes, cancer_diseases)
    count_filter(os.path.join(args.output_dir, f"Methylation_sig_promoter_out_number.txt"),
                 promoter, promoter_genes, cancer_diseases)
    count_filter(os.path.join(args.output_dir, f"Methylation_sig_body_out_number_pandisease.txt"),
                 body, body_genes, diseases)
    count_filter(os.path.join(args.output_dir, f"Methylation_sig_promoter_out_number_pandisease.txt"),
                 promoter, promoter_genes, diseases)

    # long-form labels
    label_path = os.path.join(args.output_dir, "Sig_label_out.txt")
    with open(label_path, "w") as fh:
        fh.write("geneid\tvalue\tdataset\tsig_label\tlocus\tcondition\n")
        for d in diseases:
            for region_dict, locus in ((body, "body"), (promoter, "promoter")):
                for g, (direction, ev) in region_dict[d["name"]].items():
                    if not keep(g):
                        continue
                    fh.write(f"{g}\t{ev}\t{d['dataset']}\tSig\t{locus}\t{d['name']}\n")

    # database-import pivot: one body row + one promoter row per gene
    import_path = os.path.join(args.output_dir, "Sig_label_out_database_import.txt")
    disease_names = [d["name"] for d in diseases]
    with open(import_path, "w") as fh:
        fh.write("Gene\t" + "\t".join(disease_names) + "\tLocation\n")
        for g in all_genes:
            for region_dict, locus in ((body, "body"), (promoter, "promoter")):
                cells = []
                for d in diseases:
                    hit = region_dict[d["name"]].get(g)
                    cells.append(direction_label(hit[0]) if hit else "NA")
                fh.write(g + "\t" + "\t".join(cells) + f"\t{locus}\n")

    print(f"done: {len(all_genes)} genes -> {args.output_dir}")


if __name__ == "__main__":
    main()
