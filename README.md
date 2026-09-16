# lncbookv3-processing

Nextflow DSL2 pipelines for LncBook v3 genomic annotation processing.

## Overview

This project maps external genomic resources onto LncBook v3 lncRNA annotations:

- ClinVar variants with definite clinical labels
- COSMIC pathogenic variants filtered by FATHMM-MKL score
- GWAS Catalog genome-wide significant associations
- SmProt small protein coordinates
- methylation, ceRNA and lncRNA conservation analyses

The pipeline validates inputs, normalizes records, extracts reusable GTF/GFF3/BED features, validates intervals before intersection, intersects coordinates with `bedtools`, processes pslMap-derived conservation evidence, and exports annotation tables, reports, SQL files, UCSC Genome Browser BED tracks and database-import tables.

## Project Layout

- [main.nf](main.nf): top-level Nextflow DSL2 entry point.
- [nextflow.config](nextflow.config): defaults, profiles and container settings.
- [nextflow_schema.json](nextflow_schema.json): parameter schema.
- [modules/local/](modules/local): reusable nf-core-style local modules, one process per tool.
- [modules/nf-core/](modules/nf-core): vendored nf-core modules for standard format conversions.
- [subworkflows/local/](subworkflows/local): branch-level subworkflows composed from modules.
- [tests/data/](tests/data): small smoke-test fixtures.
- [scripts/](scripts): Python implementations used by the local modules.
- [Dockerfile](Dockerfile): optional custom images for UCSC kent tools, samtools and ceRNA prediction binaries.
- [docs/nextflow_pipeline.md](docs/nextflow_pipeline.md): Nextflow-specific architecture and usage notes.
- [docs/conservation_pipeline.md](docs/conservation_pipeline.md): lncRNA conservation analysis design.

## Workflow Structure

The reference branch always runs:

1. validate the LncBook annotation,
2. extract `gene/transcript/exon/intron` BED-like features,
3. validate the extracted BED intervals.

Optional branches are controlled by `params.run_*` flags:

- `run_clinvar`: ClinVar validation, label filtering, BED validation, intersection, annotation, summary and export.
- `run_cosmic`: COSMIC validation, FATHMM-MKL filtering, BED validation, intersection, annotation, summary and export.
- `run_gwas_catalog`: GWAS Catalog validation, genome-wide-significance filtering, BED validation, intersection, annotation, summary and export.
- `run_smprot`: SmProt validation, coordinate normalization, in-transcript mapping, summary and export.
- `run_methylation`: manifest-driven methylation analysis.
- `run_cerna`: manifest-driven ceRNA integration.
- `run_conservation`: manifest-driven lncRNA conservation analysis.

## Reusable Modules

The project now follows an nf-core-inspired layout. Every tool lives in its own directory under `modules/local/<tool>/` with:

- `main.nf`: the DSL2 process,
- `meta.yml`: module metadata,
- `environment.yml`: conda dependencies for standalone reuse.

Branch logic lives in `subworkflows/local/<branch>/main.nf`. This keeps low-level tools testable and reusable while keeping the top-level pipeline readable.

BEDOPS converters are exposed as independent modules: `gtf2bed`, `gff2bed`, `vcf2bed`, `sam2bed`, `bam2bed`, `psl2bed`, `rmsk2bed` and `wig2bed`. They use `quay.io/biocontainers/bedops:2.4.42--hd6d6fdc_1` by default, or the module conda environment with `bedops=2.4.42`.

Standard format conversions and validations are provided by nf-core modules vendored under `modules/nf-core/` (source: [nf-core/modules](https://github.com/nf-core/modules) @ `0befdd9`, MIT), replacing the former custom conversion modules:

- BED to bigBed: `ucsc/bedtobigbed` (sort first with `bedtools/sort`),
- bedGraph to bigWig: `ucsc/bedgraphtobigwig`; WIG to bigWig: `ucsc/wigtobigwig`; out-of-bounds clipping: `ucsc/bedclip`,
- BED to GFF3: `agat/convertbed2gff`; GFF3 to BED12: `agat/convertgff2bed`; GFF/GTF to GFF3/GTF: `agat/convertspgxf2gxf` and `agat/convertspgff2gtf`; GTF to BED: `bedops/gtf2bed`; BAM/GFF/GTF/GVF/PSL to BED: `bedops/convert2bed`; GTF to genePred: `ucsc/gtftogenepred`; GFF/GTF validate and convert: `gffread`,
- coverage tracks: `bedtools/genomecov` (BAM/BED to bedGraph) and `deeptools/bamcoverage` (BAM to bigWig/bedGraph); BAM to BED12: `bedtools/bamtobed`,
- VCF/BCF: `bcftools/view`, `bcftools/query`; BED to VCF: `bedgovcf`,
- SAM/BAM/CRAM conversion, sorting and indexing: `samtools/view`, `samtools/sort` and `samtools/index`,
- validation: `gt/gff3validator` (GFF3), `htsnimtools/vcfcheck` (VCF), `samtools/quickcheck` (BAM/CRAM).

See [docs/nextflow_pipeline.md](docs/nextflow_pipeline.md) for import statements, recommended compositions and `ext.args` examples.

## Run

Run the reference branch only:

```bash
nextflow run . --gtf data/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf
```

Run the bundled smoke profile:

```bash
nextflow run . -profile test
```

Run with conda environments:

```bash
nextflow run . -profile test,conda
```

Run with Docker containers:

```bash
nextflow run . -profile docker --gtf data/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf
```

Enable branches by supplying both the flag and input path, for example:

```bash
nextflow run . \
  --gtf data/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf \
  --run_clinvar true \
  --clinvar_vcf data/variants/clinvar/clinvar_20260503.vcf.gz \
  --run_conservation true \
  --conservation_manifest scripts/conservation/conservation.example.json \
  --conservation_data_root data/conservation
```

## Outputs

Process outputs are copied to `results/LNCBOOKV3_PROCESSING/...` by default and are also available as named workflow emits for downstream reuse. Key logical outputs include:

- validated annotation and `gtf_features.bed`,
- branch-specific annotation TSVs,
- summary TSVs,
- SQL exports,
- UCSC BED tracks,
- validation/filtering reports,
- `methylation_outputs/`, `cerna_outputs/` and `conservation_outputs/` when those branches are enabled.

## Notes

- The legacy workflow implementation and job files have been removed.
- The former smoke-test fixtures were migrated to [tests/data/](tests/data).
- See [docs/nextflow_pipeline.md](docs/nextflow_pipeline.md) for module/subworkflow reuse examples.
