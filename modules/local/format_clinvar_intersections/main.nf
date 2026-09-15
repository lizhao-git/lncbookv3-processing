nextflow.enable.dsl = 2

process FORMAT_CLINVAR_INTERSECTIONS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path('gtf_clinvar_annotations.tsv'), emit: annotations

    script:
    """
    python3 ${projectDir}/scripts/variant_pipeline/format_intersections.py --input-intersections ${raw} --output-tsv gtf_clinvar_annotations.tsv
    """
}
