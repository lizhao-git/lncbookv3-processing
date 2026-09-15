nextflow.enable.dsl = 2

process VALIDATE_SMPROT_TSV {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(smprot_files)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated
    tuple val(meta), path(report_name), emit: report

    script:
    def inputs = smprot_files.collect { "--input-file ${it}" }.join(' ')
    """
    python3 ${projectDir}/scripts/smprot_pipeline/validate_smprot_tsv.py \\
        ${inputs} \\
        --output-tsv ${output_name} \\
        --report ${report_name}
    """
}
