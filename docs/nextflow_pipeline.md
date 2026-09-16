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

## Format conversions and validations via nf-core modules

Standard format conversions and validations use nf-core modules vendored under `modules/nf-core/` (copied from [nf-core/modules](https://github.com/nf-core/modules), MIT licence; `main.nf`, `environment.yml` and `meta.yml` only). Each module reads extra tool arguments from `task.ext.args` and the output prefix from `task.ext.prefix`. 31 modules are vendored, pinned to nf-core/modules commit [`0befdd9`](https://github.com/nf-core/modules/commit/0befdd9db2975ee04a97decc1672a4304387e9a7).

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

### Conversion and validation module library

Annotation and interval conversions:

| Module | Direction |
| --- | --- |
| `bedops/gtf2bed` | GTF → BED |
| `bedops/convert2bed` | BAM/GFF/GTF/GVF/PSL → BED (input format inferred from the file extension) |
| `agat/convertgff2bed` | GFF3 → BED12 |
| `ucsc/gtftogenepred` | GTF → genePred; with `ext.args = '-genePredExt -geneNameAsName2'` a refFlat file is produced as well |
| `gffread` | validate, filter and convert GFF/GTF; select the output with `ext.args` (`-T` → GTF, `--bed` → BED, `-w`/`-x`/`-y` → FASTA, default GFF3) |

Coverage tracks (alignment/intervals → bedGraph, bigWig, BED):

| Module | Direction |
| --- | --- |
| `ucsc/wigtobigwig` | WIG → bigWig (needs `chrom.sizes`) |
| `ucsc/bedclip` | clip bedGraph records that overrun `chrom.sizes` (typically before `ucsc/bedgraphtobigwig`) |
| `bedtools/genomecov` | BAM or BED intervals → coverage bedGraph (`-bg`), histogram or per-base reports, with optional built-in sort |
| `deeptools/bamcoverage` | BAM (+ index, optional FASTA and blacklist) → bigWig, or bedGraph with `--outFileFormat bedgraph` |
| `bedtools/bamtobed` | BAM → BED12 |
| `modkit/bedmethyltobigwig` | bedMethyl (ONT modkit) → bigWig (needs `chrom.sizes`) |

VCF/BCF conversions:

| Module | Direction |
| --- | --- |
| `bcftools/view` | subset/filter VCF/BCF by regions, targets or samples; output type via `ext.args` (`-Oz`, `-Ob`, ...) |
| `bcftools/query` | extract VCF/BCF fields into a table (`-f` format string via `ext.args`, file suffix via `ext.suffix`, default `txt`) |
| `bedgovcf` | BED + YAML config + FASTA index → bgzipped VCF |
| `gvcftools/extractvariants` | extract variants in a region list from a VCF/BCF |

Liftover, compression and indexing:

| Module | Direction |
| --- | --- |
| `ucsc/liftover` | BED + chain file → lifted and unlifted BED |
| `picard/liftovervcf` | VCF + chain + reference FASTA → lifted VCF |
| `htslib/bgziptabix` | bgzip (or decompress) a file and optionally build a tabix/CSI index |
| `tabix/bgzip` | bgzip compression with optional tabix index |

Generic format validation (complements the data-specific `validate_*` local modules used by the pipeline branches):

| Module | Checks |
| --- | --- |
| `gt/gff3validator` | strict GFF3 validation (GenomeTools); emits a `*.success.log` or `*.error.log` |
| `htsnimtools/vcfcheck` | a VCF against a background VCF (e.g. gnomAD) → TSV report |
| `samtools/quickcheck` | BAM/CRAM header and EOF integrity; the task fails on corrupt files |

Validate a GFF3 file with the bundled fixture:

```nextflow
include { GT_GFF3VALIDATOR } from './modules/nf-core/gt/gff3validator/main'

workflow {
    gff3_ch = Channel.of(tuple([id: 'sample'], file('tests/data/inputs/sample.gff3')))
    GT_GFF3VALIDATOR(gff3_ch)
}
```

### Cross-assembly and multi-species liftover

Cross-assembly liftover is available in three flavours:

- `ucsc/liftover` (vendored nf-core module): one BED file plus one chain file per task (`liftOver` under the hood, kent v482 container).
- `modules/local/liftover`: wraps `scripts/format_convert/liftover.py`. In addition to the vendored module it emits the unmapped BED, exposes `--min-match` (default `0.95`) and writes a mapping-rate report (`input_records`, `mapped_records`, `unmapped_records`, `mapping_rate`). Runs in the kent container (`container_kent`).
- `modules/local/liftover_multi`: wraps `scripts/format_convert/liftover_multi.py`, a JSON-manifest driver for many species/assemblies in one task. The manifest registers chain files per `(species, from, to)` and lists datasets referencing them; see [scripts/format_convert/liftover.example.json](../scripts/format_convert/liftover.example.json).

Manifest layout:

```json
{
  "min_match": 0.95,
  "chains": [
    {"species": "human", "from": "hg19", "to": "hg38", "chain": "chains/hg19ToHg38.over.chain.gz"},
    {"species": "mouse", "from": "mm39", "to": "mm10", "chain": "chains/mm39ToMm10.over.chain.gz"}
  ],
  "datasets": [
    {"name": "human_lncrna", "species": "human", "from": "hg19", "to": "hg38", "input": "human/lncrna.bed"},
    {"name": "mouse_other", "species": "mouse", "from": "mm39", "to": "mm10", "input": "mouse/other.bed", "chain": "chains/custom.chain.gz", "min_match": 0.9}
  ]
}
```

Datasets may override the global `min_match` and point at a dataset-specific `chain`; otherwise the chain registry is looked up by `(species, from, to)`. Relative `input`/`chain` paths are resolved against `--data-root`. Each dataset produces `<name>.lifted.bed` and `<name>.unmapped.bed` under `--output-dir` (`liftover_outputs`), and the summary `liftover_report.tsv` records `dataset, species, from, to, chain, status, input_records, mapped_records, unmapped_records, mapping_rate, message`. A failing dataset is recorded and does not abort the remaining ones, but the task exits non-zero if any dataset failed.

### Local format validators

Formats without an nf-core validation module are covered by local modules that wrap `scripts/format_convert/validate_*.py` (stdlib only; each validates and then copies the input through):

| Module | Checks |
| --- | --- |
| `validate_bigtrack` | bigWig/bigBed via kent `bigWigInfo`/`bigBedInfo` (`version`, `chromCount`, `basesCovered`, item statistics); runs in the kent container |
| `validate_bedgraph` | 4 columns, integer ranges with `start < end`, finite values, ordering warnings |
| `validate_wig` | `fixedStep`/`variableStep` state machine, implicit positions, 1-based ordinates, ordering warnings |
| `validate_genepred` | genePred vs refFlat auto-detection, exon count/list consistency, exon and CDS bounds |
| `validate_sam` | headers (@HD first, @SQ SN/LN), 11 mandatory fields, CIGAR/SEQ/QUAL consistency, optional tag syntax — no samtools needed |
| `validate_chain` | chain header fields, strand/size checks, block accounting against target/query spans, size-only terminator |
| `validate_fasta` | IUPAC alphabet (warnings for other characters), duplicate headers, empty sequences, length statistics |

`validate_chain` is also useful as a pre-flight check before `liftover`/`liftover_multi` runs, and `validate_bedgraph`/`validate_wig` before `bedGraphToBigWig`/`wigToBigWig`.

### Not available as nf-core modules

The nf-core/modules catalogue (checked at commit `0befdd9`) has no modules for these conversions and checks; keep using the kent binaries or the remaining local modules:

- bigWig → bedGraph/text and bigBed → BED: keep using the local `bigtrack_to_text` module (`bigWigToBedGraph` / `bigBedToBed`).
- Generic BED, GTF and annotation validators: the local `validate_bed`, `validate_gtf` and `validate_annotation` modules keep that role; bedGraph, WIG, genePred/refFlat, SAM, chain, FASTA and bigWig/bigBed now have their own local validators (`validate_bedgraph`, `validate_wig`, `validate_genepred`, `validate_sam`, `validate_chain`, `validate_fasta`, `validate_bigtrack`).
- BED → GATK interval_list: `gatk4/bedtointervallist`, `gatk4/intervallisttobed` and `picard/bedtointervallist` exist upstream but are GATK-specific; vendor them only if the pipeline adopts GATK tooling.

Liftover and compression are now vendored (`ucsc/liftover`, `picard/liftovervcf`, `htslib/bgziptabix`, `tabix/bgzip`). Related nf-core modules that remain unvendored because they are tool-specific variants or other storage formats: `tabix/tabix` and `samtools/bgzip` (superseded by `htslib/bgziptabix` for compress-and-index); `ea-utils/gtf2bed`, `gtfsort`, `gffcompare`; `biscuit/vcf2bed`, `svtk/vcf2bed` (tool-specific VCF → BED); `plink2/vcf2bgen`, `vcf2db`, `vcf2maf`, `vcf2zarr`, `bio2zarr/vcf2zarrconvert` (VCF → other storage/report formats).
