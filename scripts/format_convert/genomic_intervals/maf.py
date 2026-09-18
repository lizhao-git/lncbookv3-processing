"""Validate TCGA MAF (Mutation Annotation Format) files (stdlib only).

MAF is the tab-delimited variant table produced by ``vcf2maf`` and consumed
by the ``last/mafconvert`` family. Comment lines starting with ``#`` (e.g.
``#version 2.4``) are allowed before the header; the first non-comment line
must be the column header. Core coordinate/allele columns are hard errors
when missing, gene and annotation columns are softer checks because
positional-only MAFs are common.
"""
import re

from .io import open_text

CORE_COLUMNS = (
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
)
EXPECTED_COLUMNS = ("Hugo_Symbol", "Strand", "NCBI_Build", "Center", "Entrez_Gene_Id")

KNOWN_TYPES = {"SNP", "DNP", "TNP", "ONP", "INS", "DEL"}
MNP_LENGTHS = {"DNP": 2, "TNP": 3, "ONP": 4}

KNOWN_CLASSES = {
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent",
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Splice_Site", "Splice_Region", "Translation_Start_Site",
    "Start_Codon_SNP", "Start_Codon_Ins", "Start_Codon_Del",
    "Stop_Codon_SNP", "Stop_Codon_Ins", "Stop_Codon_Del",
    "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
    "RNA", "Intron", "IGR", "lincRNA", "Targeted_Region",
    "3'UTR", "5'UTR", "3'Flank", "5'Flank",
}

ALLELE_RE = re.compile(r"^[ACGTNacgtn]+$")
MISSING_ALLELES = ("", "-", ".")


def _parse_int(errors, line_no, label, token, minimum=1):
    try:
        value = int(token)
    except ValueError:
        errors.append((line_no, f"{label} is not an integer: {token!r}"))
        return None
    if value < minimum:
        errors.append((line_no, f"{label} must be >= {minimum}, got {value}"))
        return None
    return value


def validate_maf(input_path: str, output_path: str, report_path: str):
    errors = []
    warnings = []
    compression = "plain"
    records = 0
    header = None
    col_index = {}
    samples = set()
    genes = set()
    chroms = []
    variant_types = {}
    duplicates = 0
    missing_gene = 0
    chr_prefixed = 0
    bare_chroms = 0
    seen = {}
    first_data_line = None

    with open_text(input_path, preferred_exts=(".maf", ".maf.gz", ".txt")) as (fh, detected):
        compression = detected
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line.startswith("#"):
                if first_data_line is not None:
                    warnings.append((line_no, "comment line after data rows"))
                continue

            if header is None:
                header = line.split("\t")
                for name in set(header):
                    if header.count(name) > 1:
                        warnings.append((line_no, f"duplicate header column {name!r}"))
                for name in CORE_COLUMNS:
                    if name not in header:
                        errors.append((line_no, f"missing required MAF column {name!r}"))
                for name in EXPECTED_COLUMNS:
                    if name not in header:
                        warnings.append((line_no, f"expected MAF column {name!r} is absent"))
                for name in header:
                    if name not in col_index:
                        col_index[name] = header.index(name)
                continue

            if first_data_line is None:
                first_data_line = line_no

            fields = line.split("\t")
            if len(fields) != len(header):
                errors.append((line_no, f"row has {len(fields)} columns, header has {len(header)}"))
                continue

            records += 1

            def get(name):
                idx = col_index.get(name)
                return fields[idx] if idx is not None and idx < len(fields) else ""

            chrom = get("Chromosome").strip()
            if not chrom:
                if "Chromosome" in col_index:
                    errors.append((line_no, "Chromosome is empty"))
            else:
                if chrom not in chroms:
                    chroms.append(chrom)
                if chrom.lower().startswith("chr"):
                    chr_prefixed += 1
                else:
                    bare_chroms += 1

            start_i = end_i = None
            if "Start_Position" in col_index:
                start_i = _parse_int(errors, line_no, "Start_Position", get("Start_Position"))
            if "End_Position" in col_index:
                end_i = _parse_int(errors, line_no, "End_Position", get("End_Position"))
            if start_i is not None and end_i is not None and start_i > end_i:
                errors.append((line_no, f"Start_Position ({start_i}) > End_Position ({end_i})"))

            strand = get("Strand").strip()
            if strand and strand not in ("+", "-", "."):
                warnings.append((line_no, f"Strand {strand!r} is not '+', '-' or '.'"))

            variant_type = get("Variant_Type").strip()
            ref = get("Reference_Allele").strip()
            allele1 = get("Tumor_Seq_Allele1").strip()
            allele2 = get("Tumor_Seq_Allele2").strip()

            if "Variant_Type" in col_index and not variant_type:
                errors.append((line_no, "Variant_Type is empty"))
            elif variant_type:
                variant_types[variant_type] = variant_types.get(variant_type, 0) + 1
                if variant_type not in KNOWN_TYPES:
                    warnings.append((line_no, f"unrecognised Variant_Type {variant_type!r}"))
                for label, token in (("Reference_Allele", ref), ("Tumor_Seq_Allele2", allele2)):
                    if token not in MISSING_ALLELES and not ALLELE_RE.match(token):
                        warnings.append((line_no, f"{label} {token!r} is not a DNA allele"))
                if variant_type == "SNP":
                    if len(ref) != 1 or len(allele2) != 1:
                        warnings.append((line_no, "SNP alleles are not single bases"))
                elif variant_type in MNP_LENGTHS:
                    expected_len = MNP_LENGTHS[variant_type]
                    if len(ref) != expected_len or len(allele2) != expected_len:
                        warnings.append(
                            (line_no, f"{variant_type} alleles are not {expected_len} bases")
                        )
                elif variant_type == "INS":
                    if not len(allele2) > len(ref):
                        warnings.append((line_no, "INS Tumor_Seq_Allele2 is not longer than Reference_Allele"))
                elif variant_type == "DEL":
                    if not len(ref) > len(allele2):
                        warnings.append((line_no, "DEL Reference_Allele is not longer than Tumor_Seq_Allele2"))

            if ref and ref not in MISSING_ALLELES and ref == allele2:
                warnings.append((line_no, "Reference_Allele equals Tumor_Seq_Allele2"))

            classification = get("Variant_Classification").strip()
            if "Variant_Classification" in col_index and not classification:
                errors.append((line_no, "Variant_Classification is empty"))
            elif classification and classification not in KNOWN_CLASSES:
                warnings.append((line_no, f"unrecognised Variant_Classification {classification!r}"))

            barcode = get("Tumor_Sample_Barcode").strip()
            if "Tumor_Sample_Barcode" in col_index and not barcode:
                errors.append((line_no, "Tumor_Sample_Barcode is empty"))
            elif barcode:
                samples.add(barcode)
                if any(char.isspace() for char in barcode):
                    warnings.append((line_no, "Tumor_Sample_Barcode contains whitespace"))

            hugo = get("Hugo_Symbol").strip()
            if hugo:
                genes.add(hugo)
            elif "Hugo_Symbol" in col_index:
                missing_gene += 1

            key = (chrom, get("Start_Position"), get("End_Position"), barcode, ref, allele2)
            if key in seen:
                duplicates += 1
                warnings.append((line_no, f"duplicate variant for {barcode} at {chrom} "
                                          f"(first seen on line {seen[key]})"))
            else:
                seen[key] = line_no

    if header is None:
        errors.append((0, "No MAF header line found"))
    if records == 0:
        errors.append((0, "No MAF data rows found"))
    if missing_gene:
        warnings.append((0, f"{missing_gene} row(s) have an empty Hugo_Symbol"))
    if chr_prefixed and bare_chroms:
        warnings.append((0, "mixed chromosome naming ('chr'-prefixed and bare names); this breaks "
                            "annotation and liftover matching"))

    with open(report_path, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"input_compression\t{compression}\n")
        report.write(f"records\t{records}\n")
        report.write(f"header_columns\t{len(header) if header else 0}\n")
        report.write(f"unique_samples\t{len(samples)}\n")
        report.write(f"unique_genes\t{len(genes)}\n")
        report.write(f"chrom_count\t{len(chroms)}\n")
        report.write("variant_types\t" + ",".join(f"{name}:{count}" for name, count in sorted(variant_types.items())) + "\n")
        report.write(f"duplicate_rows\t{duplicates}\n")
        report.write(f"error_count\t{len(errors)}\n")
        report.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            report.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            report.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("MAF validation failed. See report for details.")

    with open_text(input_path, preferred_exts=(".maf", ".maf.gz", ".txt")) as (in_fh, _), \
            open(output_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
