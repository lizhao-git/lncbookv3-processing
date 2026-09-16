nextflow.enable.dsl = 2

process VALIDATE_BIGTRACK {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_track)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated_track
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/validate_bigtrack.py \\
        --input-track ${input_track} \\
        --output-track ${output_name} \\
        --report ${report_name} \\
        ${args}
    """
}
