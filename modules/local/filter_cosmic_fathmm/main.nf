nextflow.enable.dsl = 2

process FILTER_COSMIC_FATHMM {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(validated_tsv)

    output:
    tuple val(meta), path('cosmic_pathogenic.bed'), emit: bed
    tuple val(meta), path('cosmic_pathogenic.tsv'), emit: tsv
    tuple val(meta), path('cosmic_pathogenic_filter_report.tsv'), emit: report

    script:
    """
    python3 ${projectDir}/scripts/variant_pipeline/filter_cosmic_fathmm.py \
        --input-tsv ${validated_tsv} \
        --output-bed cosmic_pathogenic.bed \
        --output-tsv cosmic_pathogenic.tsv \
        --report cosmic_pathogenic_filter_report.tsv
    """
}
