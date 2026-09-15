nextflow.enable.dsl = 2

process FORMAT_SMPROT_INTERSECTIONS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path('gtf_smprot_annotations.tsv'), emit: annotations
    tuple val(meta), path('smprot_mapping_report.tsv'), emit: report

    script:
    """
    python3 ${projectDir}/scripts/smprot_pipeline/format_smprot_intersections.py \\
        --input-intersections ${raw} \\
        --output-tsv gtf_smprot_annotations.tsv \\
        --report smprot_mapping_report.tsv
    """
}
