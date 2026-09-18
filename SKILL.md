---
name: lncbookv3_processing
description: Project-level skill for running lncbookv3 genomic annotation pipelines with Nextflow DSL2.
license: MIT
metadata:
  skill-author: lncbookv3-processing
  skill-type: project-orchestration-skill
  primary-runtime: nextflow
---

# lncbookv3-processing

## Execution

This project uses **Nextflow DSL2** as the primary workflow engine. The top-level entry point is `main.nf`; reusable local modules live in `modules/local/`, and branch-level subworkflows live in `subworkflows/local/`.

Run the smoke profile:

```bash
nextflow run . -profile test
```

Run the reference branch against the main annotation:

```bash
nextflow run . --gtf data/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf
```

Run with conda or Docker when external binaries are needed:

```bash
nextflow run . -profile test,conda
nextflow run . -profile docker --gtf data/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf
```

## Format conversions, validations and liftover

Standard conversions and validations come from nf-core modules vendored under `modules/nf-core/` (44 modules); formats without an upstream module are covered by local modules. Standalone validators (`validate_bigtrack`, `validate_bedgraph`, `validate_wig`, `validate_genepred`, `validate_sam`, `validate_chain`, `validate_fasta`, `validate_psl`, `validate_rmsk`, `validate_chrom_sizes`) and standalone kent/Python converters (`chain_to_psl`, `psl_to_chain`, `psl_to_bed`, `psl_sort`, `psl_reps`, `psl_stats`, `chain_sort`, `chain_filter`, `chain_swap`, `chain_merge_sort`, `gene_pred_to_gtf`, `gene_pred_to_bed`, `gff3_to_genepred`, `bed_to_genepred`, `bigwig_to_wig`, `bedmethyl_to_bed`) are not wired into `main.nf`; import them individually. The BED → X direction is covered by `bed_to_gtf`, `bed_to_bed12`, `bed_to_psl` and `bed_to_bedgraph`, plus the `bed_to_bigwig` subworkflow (`subworkflows/local/bed_to_bigwig`), which chains `bed_to_bedgraph` → `ucsc/bedgraphtobigwig`. Modules that use kent binaries need the `lncbookv3-kent:latest` container (built from the `Dockerfile`, which now ships 25 kent tools); the Python-only modules run in the plain Python container. Cross-assembly liftover: use `modules/local/liftover` for a single track or `modules/local/liftover_multi` for many species/assemblies from a JSON manifest (`scripts/format_convert/liftover.example.json`).

## Variant annotation onto the annotation file (GTF / GFF3 / BED)

`subworkflows/local/clinvar` takes the ClinVar VCF and an annotation file — **GTF, GFF3 or BED** (`--annotation`, format auto-detected from the extension by `subworkflows/local/annotation_features`) — and annotates variants onto gene features. GTF/GFF3 inputs go through `VALIDATE_ANNOTATION` → `EXTRACT_ANNOTATION_FEATURES` (tokens `gene transcript exon intron cds utr` → `feature_type` values `gene`, `transcript`, `exon`, `intron`, `CDS`, `5UTR`, `3UTR`); BED inputs go through `VALIDATE_BED` → `BED_TO_FEATURES` (name column canonicalised to the feature vocabulary, unknown names become `interval` rows) → `VALIDATE_BED`. Definite ClinVar labels are kept, intersected with `bedtools intersect -wa -wb`, and exported as a per-hit table plus a `feature_type` × clinical-label summary and SQL / UCSC track files. Every overlapping feature is kept — a variant inside a CDS appears once per matching `exon` and `CDS` row; `gene`/`transcript` hits are intentionally dropped from the table by `format_intersections.py`. All step reports flow through the shared subworkflow `subworkflows/local/pipeline_reports` (`PIPELINE_REPORTS`): it aggregates them into a single `pipeline_run.log` via `modules/local/aggregate_pipeline_logs` and renders `pipeline_debug_report.md` / `.html` via `modules/local/debug_report` (overview table with OK/WARN/FAIL badges per step, errors and warnings listed first, per-step detail tables; deterministic output, stdlib-only). The four variant branches `clinvar`, `cosmic`, `gwas` and `smprot` share the identical composition — `ANNOTATION_FEATURES` (gtf/gff3/bed auto-detected) + variant validation/filtering → `INTERSECT` → format → summary → export → `PIPELINE_REPORTS` — and the manifest-driven `methylation` / `conservation` branches wrap their run report through `PIPELINE_REPORTS` too (`cerna` emits no report file and is not wrapped). Their outputs are `cosmic_annotations.tsv/.sql` + `cosmic_ucsc.bed` + `cosmic_summary.tsv`, `gwas_annotations.tsv/.sql` + `gwas_ucsc.bed` + `gwas_summary.tsv` and `smprot_annotations.tsv/.sql` + `smprot_ucsc.bed` + `smprot_mapping_report.tsv` + `smprot_summary.tsv` (the old `gtf_` prefix is gone everywhere). `REFERENCE` now only validates the canonical GTF and publishes its feature intervals; no branch consumes them.

Validate a PSL against a `chrom.sizes` file, which also makes kent `pslCheck` verify target names, sizes and coordinates:

```nextflow
include { VALIDATE_PSL } from './modules/local/validate_psl/main'

workflow {
    psl_ch = Channel.of(tuple([id: 'sample'], file('tests/data/inputs/sample.psl'),
                              file('tests/data/inputs/chrom.sizes')))
    VALIDATE_PSL(psl_ch, 'validated.psl', 'psl_report.tsv')
}
```

The optional sizes file is the third element of the tuple; pass `[]` in its place to run the structural pass alone.

## Layout

- `main.nf`: top-level workflow.
- `nextflow.config`: profiles and parameters.
- `modules/local/<tool>/main.nf`: one reusable process per tool.
- `modules/nf-core/<tool>/<subtool>/main.nf`: vendored nf-core modules for standard format conversions and validations.
- `subworkflows/local/<branch>/main.nf`: branch pipelines composed from modules.
- `tests/data`: smoke-test fixtures.
- `scripts`: Python implementations used by modules.
