nextflow.enable.dsl = 2

process VALIDATE_FASTA {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(fasta)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated_fasta
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/validate_fasta.py \\
        --input-fasta ${fasta} \\
        --output-fasta ${output_name} \\
        --report ${report_name} \\
        ${args}
    """
}
