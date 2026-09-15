nextflow.enable.dsl = 2

process VALIDATE_BED {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(bed)
    val min_columns
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated_bed
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/validate_bed.py \\
        --input-bed ${bed} \\
        --output-bed ${output_name} \\
        --report ${report_name} \\
        --min-columns ${min_columns} \\
        ${args}
    """
}
