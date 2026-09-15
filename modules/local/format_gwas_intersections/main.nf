nextflow.enable.dsl = 2

process FORMAT_GWAS_INTERSECTIONS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path('gtf_gwas_annotations.tsv'), emit: annotations

    script:
    """
    python3 ${projectDir}/scripts/variant_pipeline/format_gwas_catalog_intersections.py --input-intersections ${raw} --output-tsv gtf_gwas_annotations.tsv
    """
}
