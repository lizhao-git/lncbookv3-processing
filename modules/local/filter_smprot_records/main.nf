nextflow.enable.dsl = 2

process FILTER_SMPROT_RECORDS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(validated_tsv)

    output:
    tuple val(meta), path('smprot_filtered.bed'), emit: bed
    tuple val(meta), path('smprot_filtered.tsv'), emit: tsv
    tuple val(meta), path('smprot_filtered_filter_report.tsv'), emit: report

    script:
    """
    python3 ${projectDir}/scripts/smprot_pipeline/filter_smprot_records.py \
        --input-tsv ${validated_tsv} \
        --output-bed smprot_filtered.bed \
        --output-tsv smprot_filtered.tsv \
        --report smprot_filtered_filter_report.tsv
    """
}
