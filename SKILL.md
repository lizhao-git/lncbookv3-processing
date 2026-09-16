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

Standard conversions and validations come from nf-core modules vendored under `modules/nf-core/`; formats without an upstream module are covered by local modules (`validate_bigtrack`, `validate_bedgraph`, `validate_wig`, `validate_genepred`, `validate_sam`, `validate_chain`, `validate_fasta`). Cross-assembly liftover: use `modules/local/liftover` for a single track or `modules/local/liftover_multi` for many species/assemblies from a JSON manifest (`scripts/format_convert/liftover.example.json`).

## Layout

- `main.nf`: top-level workflow.
- `nextflow.config`: profiles and parameters.
- `modules/local/<tool>/main.nf`: one reusable process per tool.
- `modules/nf-core/<tool>/<subtool>/main.nf`: vendored nf-core modules for standard format conversions and validations.
- `subworkflows/local/<branch>/main.nf`: branch pipelines composed from modules.
- `tests/data`: smoke-test fixtures.
- `scripts`: Python implementations used by modules.
