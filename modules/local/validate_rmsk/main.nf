nextflow.enable.dsl = 2

process VALIDATE_RMSK {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(rmsk)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated_rmsk
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/validate_rmsk.py \\
        --input-rmsk ${rmsk} \\
        --output-rmsk ${output_name} \\
        --report ${report_name} \\
        ${args}
    """
}
