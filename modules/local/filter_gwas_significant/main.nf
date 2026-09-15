nextflow.enable.dsl = 2

process FILTER_GWAS_SIGNIFICANT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(validated_tsv)

    output:
    tuple val(meta), path('gwas_catalog_significant.bed'), emit: bed
    tuple val(meta), path('gwas_catalog_significant.tsv'), emit: tsv
    tuple val(meta), path('gwas_catalog_significant_filter_report.tsv'), emit: report

    script:
    """
    python3 ${projectDir}/scripts/variant_pipeline/filter_gwas_catalog_significant.py \
        --input-tsv ${validated_tsv} \
        --output-bed gwas_catalog_significant.bed \
        --output-tsv gwas_catalog_significant.tsv \
        --report gwas_catalog_significant_filter_report.tsv
    """
}
