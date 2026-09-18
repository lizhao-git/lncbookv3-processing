"""Validate PLINK 1.x .map/.ped pairs (stdlib only).

The ``.map`` file holds 3-4 columns (chromosome, SNP id, optional cM and bp)
and the ``.ped`` file holds 6 mandatory columns plus two allele tokens per
variant, so the two files are checked together: the genotype token count per
sample must be ``6 + 2 x variants``.
"""
import re

from .io import open_text

CHROM_RE = re.compile(r"^(chr)?([0-9]+|[xymt]+|xy)$", re.IGNORECASE)
VALID_ALLELES = set("ACGT0IDacgtid")
KNOWN_SEX = {"0", "1", "2"}


def validate_plink(map_path: str, ped_path: str, output_map_path: str,
                   output_ped_path: str, report_path: str):
    errors = []
    warnings = []
    map_compression = "plain"
    ped_compression = "plain"

    # ---- .map ------------------------------------------------------------
    variants = 0
    dup_snps = 0
    map_chroms = []
    three_column = 0
    unsorted_map = 0
    seen_snps = {}
    last_map_chrom = None
    last_bp = None

    with open_text(map_path, preferred_exts=(".map",)) as (fh, detected):
        map_compression = detected
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            tokens = line.split()
            if len(tokens) == 3:
                three_column += 1
                chrom, snp, bp = tokens
                cM = None
            elif len(tokens) == 4:
                chrom, snp, cM, bp = tokens
            else:
                errors.append((line_no, f".map needs 4 columns (3 allowed without cM), got {len(tokens)}"))
                continue

            variants += 1

            if not chrom:
                errors.append((line_no, "chromosome is empty"))
            else:
                if chrom not in map_chroms:
                    map_chroms.append(chrom)
                if not CHROM_RE.match(chrom):
                    warnings.append((line_no, f"chromosome code {chrom!r} is not a standard PLINK code"))

            if not snp:
                errors.append((line_no, "SNP id is empty"))
            elif snp in seen_snps:
                dup_snps += 1
                warnings.append((line_no, f"duplicate SNP id {snp!r} (first seen on line {seen_snps[snp]})"))
            else:
                seen_snps[snp] = line_no

            if cM is not None:
                try:
                    cM_f = float(cM)
                    if cM_f != cM_f or cM_f < 0:
                        warnings.append((line_no, f"cM value {cM!r} is negative or not finite"))
                except ValueError:
                    errors.append((line_no, f"cM is not numeric: {cM!r}"))

            bp_i = None
            try:
                bp_i = int(bp)
            except ValueError:
                errors.append((line_no, f"bp position is not an integer: {bp!r}"))
            else:
                if bp_i < 0:
                    errors.append((line_no, f"bp must be >= 0, got {bp_i}"))
                elif bp_i == 0:
                    warnings.append((line_no, "bp position is 0"))
                elif last_bp is not None:
                    if chrom == last_map_chrom and bp_i < last_bp:
                        unsorted_map += 1
                    elif chrom != last_map_chrom and chrom in map_chroms[:-1]:
                        unsorted_map += 1
                if bp_i >= 0:
                    last_map_chrom, last_bp = chrom, bp_i

    if variants == 0:
        errors.append((0, "No variants found in the .map file"))
    if three_column:
        warnings.append((0, f"{three_column} .map row(s) use the 3-column layout (cM omitted)"))
    if unsorted_map:
        warnings.append((0, f"{unsorted_map} .map record(s) break chromosome/position ordering"))

    # ---- .ped ------------------------------------------------------------
    expected_tokens = 6 + 2 * variants
    samples = 0
    seen_people = {}
    dup_samples = 0
    odd_sex = 0
    odd_phenotype = 0
    odd_alleles = 0
    bad_genotypes = 0

    with open_text(ped_path, preferred_exts=(".ped",)) as (fh, detected):
        ped_compression = detected
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            tokens = line.split()
            if len(tokens) != expected_tokens:
                errors.append(
                    (line_no, f".ped has {len(tokens)} tokens, expected {expected_tokens} "
                              f"(6 columns + 2 x {variants} variants from the .map)")
                )
                continue

            samples += 1
            fid, iid = tokens[0], tokens[1]
            if fid == "0" or iid == "0":
                warnings.append((line_no, f"sample ID uses the missing code '0' (FID={fid}, IID={iid})"))

            person = f"{fid} {iid}"
            if person in seen_people:
                dup_samples += 1
                warnings.append((line_no, f"duplicate sample {person!r} (first seen on line {seen_people[person]})"))
            else:
                seen_people[person] = line_no

            if tokens[4] not in KNOWN_SEX:
                odd_sex += 1

            try:
                float(tokens[5])
            except ValueError:
                errors.append((line_no, f"phenotype is not numeric: {tokens[5]!r}"))
            else:
                if tokens[5] not in {"0", "1", "2", "-9"}:
                    odd_phenotype += 1

            genotypes = tokens[6:]
            for allele_a, allele_b in zip(genotypes[0::2], genotypes[1::2]):
                if len(allele_a) != 1 or len(allele_b) != 1:
                    bad_genotypes += 1
                    errors.append(
                        (line_no, f"genotype alleles must be single characters, got {allele_a!r}/{allele_b!r}")
                    )
                elif allele_a not in VALID_ALLELES or allele_b not in VALID_ALLELES:
                    odd_alleles += 1

    if samples == 0:
        errors.append((0, "No samples found in the .ped file"))
    if odd_sex:
        warnings.append((0, f"{odd_sex} sample(s) use a sex code outside 0/1/2"))
    if odd_phenotype:
        warnings.append((0, f"{odd_phenotype} sample(s) use a phenotype outside 0/1/2/-9"))
    if odd_alleles:
        warnings.append((0, f"{odd_alleles} genotype(s) use allele codes outside A/C/G/T/0/I/D"))

    with open(report_path, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"map_compression\t{map_compression}\n")
        report.write(f"ped_compression\t{ped_compression}\n")
        report.write(f"variants\t{variants}\n")
        report.write(f"samples\t{samples}\n")
        report.write(f"map_columns\t{'3' if three_column else '4'}\n")
        report.write(f"duplicate_snp_ids\t{dup_snps}\n")
        report.write(f"duplicate_samples\t{dup_samples}\n")
        report.write(f"unusual_sex_codes\t{odd_sex}\n")
        report.write(f"unusual_phenotypes\t{odd_phenotype}\n")
        report.write(f"unusual_alleles\t{odd_alleles}\n")
        report.write(f"bad_genotypes\t{bad_genotypes}\n")
        report.write(f"unsorted_map_records\t{unsorted_map}\n")
        report.write(f"error_count\t{len(errors)}\n")
        report.write(f"warning_count\t{len(warnings)}\n")
        if errors:
            report.write("errors\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in errors[:100]) + "\n")
        if warnings:
            report.write("warnings\t" + " | ".join(f"line {ln}: {msg}" for ln, msg in warnings[:100]) + "\n")

    if errors:
        raise SystemExit("PLINK validation failed. See report for details.")

    with open_text(map_path, preferred_exts=(".map",)) as (in_fh, _), \
            open(output_map_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)

    with open_text(ped_path, preferred_exts=(".ped",)) as (in_fh, _), \
            open(output_ped_path, "w", encoding="utf-8") as out_fh:
        for raw in in_fh:
            out_fh.write(raw)
