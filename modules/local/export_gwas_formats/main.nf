nextflow.enable.dsl = 2

process EXPORT_GWAS_FORMATS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(annotations)

    output:
    tuple val(meta), path('gtf_gwas_annotations.sql'), emit: sql
    tuple val(meta), path('gtf_gwas_ucsc.bed'), emit: ucsc

    script:
    """
    python3 ${projectDir}/scripts/variant_pipeline/export_gwas_catalog_annotation_formats.py \
        --input-tsv ${annotations} \
        --output-sql gtf_gwas_annotations.sql \
        --output-ucsc-bed gtf_gwas_ucsc.bed
    """
}
