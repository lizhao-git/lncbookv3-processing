nextflow.enable.dsl = 2

process VALIDATE_BEDGRAPH {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(bedgraph)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated_bedgraph
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/validate_bedgraph.py \\
        --input-bedgraph ${bedgraph} \\
        --output-bedgraph ${output_name} \\
        --report ${report_name} \\
        ${args}
    """
}
