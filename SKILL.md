---
name: lncbookv3_processing
description: Project-level skill for running lncbookv3 genomic annotation pipelines with CWL + Docker.
license: MIT
metadata:
  skill-author: lncbookv3-processing
  skill-type: project-orchestration-skill
  primary-runtime: cwl
---

# lncbookv3-processing

## Execution

This project uses **CWL** as the primary workflow engine with a single Docker runtime image. Build the image once, then run the workflow with `cwltool`:

```bash
docker build -t lncbookv3-processing:latest .
cwltool --outdir results/cwl cwl/lncbookv3.cwl cwl/lncbookv3-job.yml
```

The image bundles Python 3.11, bedtools, BEDOPS, samtools and the UCSC kent tools; Python processes and interval operations run inside it. The gate `ExpressionTool`s use JavaScript, so a `node` binary must be on `PATH`, or the Docker daemon must be running (cwltool then uses `node:alpine`).

## Examples

```bash
# Default run (ClinVar + GWAS + SmProt; COSMIC disabled)
cwltool --outdir results/cwl cwl/lncbookv3.cwl cwl/lncbookv3-job.yml

# Enable COSMIC after preparing an extracted TSV (edit the job file):
#   run_cosmic: true
#   cosmic_tsv: {class: File, path: data/variants/cosmic/cosmic_variants.tsv}

# Run only ClinVar (set the other run_* flags to false, files to null / [])
```
