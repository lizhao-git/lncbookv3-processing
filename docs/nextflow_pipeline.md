# lncbookv3-processing Nextflow pipeline

This project now uses Nextflow DSL2 with an nf-core-inspired layout:

- `main.nf`: top-level entry point.
- `nextflow.config`: profiles and default parameters.
- `modules/local/<tool>/main.nf`: one reusable process per local tool.
- `modules/local/<tool>/meta.yml`: short module metadata.
- `modules/local/<tool>/environment.yml`: conda environment for the process.
- `subworkflows/local/<branch>/main.nf`: reusable branch workflows composed from modules.
- `tests/data`: small smoke-test fixtures migrated from the previous workflow layout.

## Branches

The reference branch always runs and validates the annotation before extracting reusable BED features. Optional branches are controlled with `params.run_*` flags:

- `run_clinvar`: ClinVar VCF validation, label filtering, BED validation, intersection, annotation, summary and export.
- `run_cosmic`: COSMIC TSV validation, FATHMM-MKL filtering, intersection, annotation, summary and export.
- `run_gwas_catalog`: GWAS Catalog TSV validation, significance filtering, intersection, annotation, summary and export.
- `run_smprot`: SmProt coordinate validation, filtering, in-transcript mapping, summary and export.
- `run_methylation`: JSON-manifest-driven methylation analysis.
- `run_cerna`: JSON-manifest-driven ceRNA integration.
- `run_conservation`: JSON-manifest-driven lncRNA conservation analysis.

## Examples

Run the reference branch only:

```bash
nextflow run . --gtf data/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf
```

Run the bundled smoke profile:

```bash
nextflow run . -profile test,conda
```

Run with Docker containers where available:

```bash
nextflow run . -profile docker --gtf data/LncBook_v3_hg38.lncRNAs_attr_normalized.gtf
```

## Single-module reuse

Every module can be imported from another DSL2 workflow, for example:

```nextflow
include { VALIDATE_BED } from './modules/local/validate_bed/main'

workflow {
    bed_ch = Channel.of(tuple([id: 'example'], file('tests/data/inputs/sample.bed')))
    VALIDATE_BED(bed_ch, 3, 'validated.bed', 'validate_bed_report.tsv')
}
```

BEDOPS converters are available as local modules (`gtf2bed`, `gff2bed`, `vcf2bed`, `sam2bed`, `bam2bed`, `psl2bed`, `rmsk2bed`, `wig2bed`) and use the BioContainers BEDOPS image or the module conda environment.

## Format conversions via nf-core modules

Standard format conversions use nf-core modules vendored under `modules/nf-core/` (copied from [nf-core/modules](https://github.com/nf-core/modules), MIT licence; `main.nf`, `environment.yml` and `meta.yml` only). Each module reads extra tool arguments from `task.ext.args` and the output prefix from `task.ext.prefix`.

The packed converter modules (`bed_to_bigbed` and `bedgraph_to_bigwig`) are gone: `bedToBigBed` and `bedGraphToBigWig` expect records grouped per chromosome with ascending starts (`bedGraphToBigWig` offers no `-sort` flag), so sort first with `bedtools/sort` and pass the same `chrom.sizes` to the converter:

```nextflow
include { BEDTOOLS_SORT }          from './modules/nf-core/bedtools/sort/main'
include { UCSC_BEDTOBIGBED }       from './modules/nf-core/ucsc/bedtobigbed/main'
include { UCSC_BEDGRAPHTOBIGWIG }  from './modules/nf-core/ucsc/bedgraphtobigwig/main'
include { AGAT_CONVERTBED2GFF }    from './modules/nf-core/agat/convertbed2gff/main'

workflow {
    bed_ch      = Channel.of(tuple([id: 'sample'], file('tests/data/inputs/sample.bed')))
    bedgraph_ch = Channel.of(tuple([id: 'sample'], file('tests/data/inputs/sample.bedgraph')))
    chrom_sizes = file('tests/data/inputs/chrom.sizes')

    BEDTOOLS_SORT(bed_ch, chrom_sizes)
    UCSC_BEDTOBIGBED(BEDTOOLS_SORT.out.sorted, chrom_sizes, [])

    BEDTOOLS_SORT(bedgraph_ch, chrom_sizes)
    UCSC_BEDGRAPHTOBIGWIG(BEDTOOLS_SORT.out.sorted, chrom_sizes)

    AGAT_CONVERTBED2GFF(bed_ch)
}
```

Set module arguments in the pipeline config (or per include site):

```nextflow
process {
    withName: 'UCSC_BEDTOBIGBED'    { ext.args = '-type=bed6' }  // optional: standard BED3-BED12 layouts are auto-detected; -type is only needed for non-standard bedPlus layouts
    withName: 'AGAT_CONVERTBED2GFF' { ext.args = '--source lncbookv3 --primary_tag lncRNA' }
}
```

Tool behaviour verified against the kent v482 BioContainers images shipped with the modules:

- `bedToBigBed` auto-detects standard BED3-BED12 field counts; its own `-sort` flag also works, and the chromosome order in the input does not have to match `chrom.sizes`. Pre-sorting with `BEDTOOLS_SORT` is still preferred so the sorted channel can be reused by other steps.
- `bedGraphToBigWig` in v482 has neither a `-clip` nor a `-sort` option: coordinates that overrun `chrom.sizes` fail the task, so clip upstream if needed.

The former custom conversion modules were replaced as follows:

| Former local module | nf-core modules now used |
| --- | --- |
| `bed_to_bigbed` (sorting, BED type inference, report) | `bedtools/sort` then `ucsc/bedtobigbed` |
| `bedgraph_to_bigwig` (sorting, `--clip`, report) | `bedtools/sort` then `ucsc/bedgraphtobigwig` (kent v482 has no `-clip`; clip upstream instead) |
| `bed_to_gff3` | `agat/convertbed2gff` |
| `convert_annotation_format` | `agat/convertspgxf2gxf` (to GFF3) / `agat/convertspgff2gtf` (to GTF) |
| `convert_alignment_format` | `samtools/view` (plus `samtools/sort`, `samtools/index`) |
