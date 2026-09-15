nextflow.enable.dsl = 2

process VALIDATE_GTF {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(gtf)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated_gtf
    tuple val(meta), path(report_name), emit: report

    script:
    """
    python3 ${projectDir}/scripts/format_convert/validate_gtf.py \\
        --input-gtf ${gtf} \\
        --output-gtf ${output_name} \\
        --report ${report_name}
    """
}
