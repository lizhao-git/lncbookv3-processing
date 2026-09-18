nextflow.enable.dsl = 2

process FORMAT_COSMIC_INTERSECTIONS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(raw)
    val output_prefix

    output:
    tuple val(meta), path("${output_prefix}_annotations.tsv"), emit: annotations

    script:
    """
    python3 ${projectDir}/scripts/variant_pipeline/format_cosmic_intersections.py --input-intersections ${raw} --output-tsv ${output_prefix}_annotations.tsv
    """
}
