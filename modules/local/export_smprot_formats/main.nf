nextflow.enable.dsl = 2

process EXPORT_SMPROT_FORMATS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(annotations)
    val output_prefix

    output:
    tuple val(meta), path("${output_prefix}_annotations.sql"), emit: sql
    tuple val(meta), path("${output_prefix}_ucsc.bed"), emit: ucsc

    script:
    """
    python3 ${projectDir}/scripts/smprot_pipeline/export_smprot_annotation_formats.py \
        --input-tsv ${annotations} \
        --output-sql ${output_prefix}_annotations.sql \
        --output-ucsc-bed ${output_prefix}_ucsc.bed
    """
}
