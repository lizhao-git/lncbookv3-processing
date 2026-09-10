# lncbookv3-processing

CWL + Docker pipelines for LncBook v3 genomic annotation processing.

## Overview

This project maps external genomic resources onto LncBook v3 lncRNA annotations:

- ClinVar variants with definite clinical labels
- COSMIC pathogenic variants filtered by FATHMM-MKL score
- GWAS Catalog genome-wide significant associations
- SmProt small protein coordinates

The pipeline validates inputs, normalizes records, extracts lncRNA GTF/GFF3 features, validates BED intervals before intersection, intersects coordinates with `bedtools`, and exports annotation tables, summary reports, SQL files, and UCSC Genome Browser BED tracks.

The primary workflow engine is **CWL** (Common Workflow Language) with a single **Docker** runtime image.

## Project Layout

- [cwl/lncbookv3.cwl](cwl/lncbookv3.cwl): top-level CWL workflow.
- [cwl/lncbookv3-job.yml](cwl/lncbookv3-job.yml): example job input.
- [cwl/workflows/](cwl/workflows): per-branch sub-workflows (reference/clinvar/cosmic/gwas/smprot/methylation/cerna/cerna_predict).
- [cwl/tools/](cwl/tools): `CommandLineTool` and gate (`ExpressionTool`) definitions.
- [Dockerfile](Dockerfile): runtime image (Python 3.11 + bedtools + Bedops + samtools + UCSC kent tools + bundled scripts).
- [scripts/](scripts): Python tools used by the pipeline.
- [scripts/format_convert/](scripts/format_convert): GTF/GFF3/BED/SAM/BAM/VCF format conversion, VCF processing, UCSC track conversion, feature extraction, and validation entry points, plus the shared `genomic_intervals` parsing library.
- [scripts/variant_pipeline/](scripts/variant_pipeline): ClinVar/COSMIC/GWAS Catalog variant-processing scripts (validation, filtering, intersection formatting, summarization, export).
- [scripts/smprot_pipeline/](scripts/smprot_pipeline): SmProt small-protein coordinate validation, filtering, intersection formatting, summarization, and export scripts.
- [scripts/methylation/](scripts/methylation): refactored methylation analysis pipeline (Python 3).
- [scripts/cerna_pipeline/](scripts/cerna_pipeline): refactored ceRNA analysis pipeline (Python 3).
- [data/](data): local input data.

## Reusable Genomic Interval Layer

Annotation and interval handling is factored out from the reference branch:

- `scripts/format_convert/genomic_intervals/annotation.py`: parses GTF/GFF3 attributes, validates 9-column annotations, copies validated annotation text, and extracts `gene/transcript/exon/intron` intervals.
- `scripts/format_convert/genomic_intervals/bed.py`: validates BED coordinate files as 0-based half-open intervals.
- `scripts/format_convert/validate_annotation.py`, `scripts/format_convert/extract_annotation_features.py`, and `scripts/format_convert/validate_bed.py`: command-line entry points used by CWL.
- `scripts/format_convert/validate_gtf.py` and `scripts/format_convert/extract_gtf_features.py`: compatibility wrappers for existing GTF-only calls.

Reusable CWL modules:

- [cwl/workflows/annotation_intervals.cwl](cwl/workflows/annotation_intervals.cwl): validate GTF/GFF3 and emit BED-like feature intervals.
- [cwl/workflows/bed_validation.cwl](cwl/workflows/bed_validation.cwl): validate existing BED files before intersection or downstream annotation.
- [cwl/tools/validate_annotation.cwl](cwl/tools/validate_annotation.cwl), [cwl/tools/extract_annotation_features.cwl](cwl/tools/extract_annotation_features.cwl), and [cwl/tools/validate_bed.cwl](cwl/tools/validate_bed.cwl): low-level reusable tools.

## Format Conversion, VCF and Track Utilities

Standalone tools in [scripts/format_convert/](scripts/format_convert), each with a 1:1 `CommandLineTool` in [cwl/tools/](cwl/tools). They are not wired into the main workflow; run them individually.

- Annotation conversion: `convert_annotation_format.py` (GTF <-> GFF3), `bed_to_gff3.py` (BED -> GFF3).
- Variant conversion: `vcf_to_bed.py` (VCF -> BED).
- Alignment conversion: `convert_alignment_format.py` (SAM <-> BAM via samtools, optional sort/index), `alignment_to_bed.py` (SAM/BAM -> BED6, optional per-block rows and MAPQ filter).
- VCF processing: `validate_vcf.py` (header/record/INFO/FORMAT validation with copy-through), `vcf_filter.py` (region/QUAL/DP/allele-type/FILTER), `vcf_sort.py` (contig + position), `vcf_stats.py` (record/allele counts, Ts/Tv), `vcf_split.py` (per-contig / per-sample), `vcf_normalize.py` (split multiallelics, subset Number=A/R/G values).
- UCSC tracks: `make_ucsc_track.py` (add/replace a track line), `bed_to_bigbed.py` and `bedgraph_to_bigwig.py` (via kent tools, rows sorted by `chrom.sizes`), `bigtrack_to_text.py` (bigBed -> BED, bigWig -> bedGraph).

Example (uses the runtime image, like the workflows):

```bash
cwltool cwl/tools/vcf_stats.cwl --vcf data/variants/clinvar/clinvar_20260503.vcf.gz
```

## Inputs

Declared by the top-level workflow ([cwl/lncbookv3.cwl](cwl/lncbookv3.cwl)):

| Input | Type | Description |
| --- | --- | --- |
| `gtf` | `File` | LncBook v3 GTF. |
| `run_clinvar` / `clinvar_vcf` | `boolean` / `File?` | Enable the ClinVar branch and supply its VCF. |
| `run_cosmic` / `cosmic_tsv` | `boolean` / `File?` | Enable the COSMIC branch and supply an extracted TSV (disabled by default). |
| `run_gwas_catalog` / `gwas_tsv` | `boolean` / `File?` | Enable the GWAS Catalog branch and supply its TSV. |
| `run_smprot` / `smprot_tsv` | `boolean` / `File[]` | Enable the SmProt branch and supply one or more coordinate files. |
| `run_methylation` / `methylation_manifest` / `methylation_data_root` | `boolean` / `File?` / `Directory?` | Enable the methylation branch; supply a JSON manifest and the data-root Directory that anchors its relative paths. |
| `run_cerna` / `cerna_manifest` / `cerna_data_root` | `boolean` / `File?` / `Directory?` | Enable the ceRNA branch; supply a JSON manifest and the data-root Directory that anchors its relative paths. |

Set a branch's file to `null` (or `[]` for `smprot_tsv`) and its `run_*` flag to `false` to skip it. See [cwl/lncbookv3-job.yml](cwl/lncbookv3-job.yml).

## Workflows

ClinVar branch:

1. Validate GTF/GFF3 annotation.
2. Validate ClinVar VCF.
3. Extract and validate `gene/transcript/exon/intron` BED features from annotation.
4. Keep ClinVar variants with one of: `Benign`, `Pathogenic`, `Affects`, `Drug response`, `Protective`, `Risk factor`.
5. Validate filtered BED intervals and intersect annotation features with filtered ClinVar variants.
6. Export annotation TSV, summary TSV, SQL, and UCSC BED.

COSMIC branch:

1. Validate COSMIC TSV.
2. Keep variants with `FATHMM-MKL > 0.7`.
3. Validate filtered BED intervals and intersect with annotation features.
4. Export annotation TSV, summary TSV, SQL, and UCSC BED.

GWAS Catalog branch:

1. Validate GWAS Catalog TSV.
2. Keep associations with `p-value < 5e-8`.
3. Validate filtered BED intervals and intersect with annotation features.
4. Export annotation TSV, summary TSV, SQL, and UCSC BED.

SmProt branch:

1. Validate and normalize one or more SmProt coordinate files.
2. Convert valid records to BED.
3. Validate filtered BED intervals and intersect with annotation features.
4. Keep mappings entirely contained in one transcript context and uniquely mapped by `protein_id`.
5. Export annotation TSV, mapping report, summary TSV, SQL, and UCSC BED.

Methylation branch:

1. Normalize per-dataset methylation data to 0-based BED (bigwig/bigbed/bismark/bed), with optional hg19→hg38 liftover.
2. Compute per-gene mean methylation over gene body and promoter regions.
3. Differential methylation testing (Wilcoxon + FDR, or fold-change consistency).
4. Aggregate cross-dataset significance into a gene × disease matrix, labels and a database-import table.

See [scripts/methylation/README.md](scripts/methylation/README.md).

ceRNA branch:

0. (Optional upstream) Predict miRNA targets from mature miRNA and GTF RNA sequences with miRanda / RNAhybrid / TargetScan — see [cwl/workflows/cerna_predict.cwl](cwl/workflows/cerna_predict.cwl).
1. Normalize miRanda / RNAhybrid / TargetScan interaction predictions into a unified table.
2. Intersect the three tools and refine binding sites (TargetScan site within miRanda/RNAhybrid site, plus a 6-nt sliding window).
3. Map experimentally validated ceRNA (LncRNAWiki + HGNC) and produce the final interaction table.
4. Annotate diseases via HMDD.

See [scripts/cerna_pipeline/README.md](scripts/cerna_pipeline/README.md).

## Run

### 1. Build the runtime image

```bash
docker build -t lncbookv3-processing:latest .
```

### 2. Run the workflow

Install `cwltool` (`pip install cwltool`). The gate expressions use JavaScript, so either a `node` binary must be on `PATH`, or the Docker daemon must be running (cwltool then uses a `node:alpine` container).

```bash
cwltool --outdir results/cwl cwl/lncbookv3.cwl cwl/lncbookv3-job.yml
```

Customise the run by editing [cwl/lncbookv3-job.yml](cwl/lncbookv3-job.yml) or supplying an inline job:

```bash
cwltool --outdir results/cwl cwl/lncbookv3.cwl \
  <(echo 'gtf: {class: File, path: data/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf}
run_clinvar: true
clinvar_vcf: {class: File, path: data/variants/clinvar/clinvar_20260503.vcf.gz}
run_cosmic: false
cosmic_tsv: null
run_gwas_catalog: false
gwas_tsv: null
run_smprot: false
smprot_tsv: []')
```

Run only ClinVar by setting the other three `run_*` flags to `false` (and their files to `null` / `[]`).

Run only the methylation workflow:

```bash
cwltool --outdir results/cwl/methylation \
  cwl/workflows/methylation.cwl \
  cwl/methylation-job.yml
```

Run only the ceRNA workflow:

```bash
cwltool --outdir results/cwl/cerna \
  cwl/workflows/cerna.cwl \
  cwl/cerna-job.yml
```

Run the ceRNA target-prediction stage (mature miRNA + GTF -> prediction files):

```bash
cwltool --outdir results/cwl/cerna_predict \
  cwl/workflows/cerna_predict.cwl \
  cwl/cerna_predict-job.yml
```

## Outputs

All branch results are collected under the `--outdir` directory (default `results/cwl`). Key files:

- `gtf_clinvar_site_annotations.tsv`
- `gtf_cosmic_site_annotations.tsv`
- `gtf_gwas_catalog_site_annotations.tsv`
- `gtf_smprot_site_annotations.tsv`
- `*_summary.tsv`
- `*_site_annotations.sql`
- `*_ucsc_track.bed`
- validation and filtering reports (e.g. `*_validation_report.tsv`, `*_filter_report.tsv`, `smprot_mapping_report.tsv`)
- `gtf_features_validation_report.tsv` for the reusable annotation feature BED.
- `methylation_outputs/` and `methylation_report.tsv` when the methylation branch is enabled.
- `cerna_outputs/` when the ceRNA branch is enabled (predicted interactions, experimental mapping, disease annotation).
- `cerna_predict/` when the ceRNA prediction stage is run (prepared sequences + normalized miRanda/RNAhybrid/TargetScan predictions).

Disabled branches produce empty output lists.
