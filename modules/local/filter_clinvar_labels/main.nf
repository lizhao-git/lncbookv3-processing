nextflow.enable.dsl = 2

process FILTER_CLINVAR_LABELS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(validated_vcf)

    output:
    tuple val(meta), path('clinvar_definite_labels.bed'), emit: bed
    tuple val(meta), path('clinvar_definite_labels.vcf'), emit: vcf
    tuple val(meta), path('clinvar_filter_report.tsv'), emit: report

    script:
    """
    python3 ${projectDir}/scripts/variant_pipeline/filter_clinvar_labels.py \\
        --input-vcf ${validated_vcf} \\
        --output-bed clinvar_definite_labels.bed \\
        --output-vcf clinvar_definite_labels.vcf \\
        --report clinvar_filter_report.tsv
    """
}
