nextflow.enable.dsl = 2

process CONVERT_ANNOTATION_FORMAT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(annotation)
    val from_format
    val to_format
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: annotation
    tuple val(meta), path(report_name), optional: true, emit: report

    script:
    def report_arg = report_name ? "--report ${report_name}" : ''
    """
    python3 ${projectDir}/scripts/format_convert/convert_annotation_format.py \\
        --input-annotation ${annotation} \\
        --output-annotation ${output_name} \\
        --from ${from_format} \\
        --to ${to_format} \\
        ${report_arg}
    """
}
