nextflow.enable.dsl = 2

process VALIDATE_PLINK {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(map), path(ped)
    val output_map_name
    val output_ped_name
    val report_name

    output:
    tuple val(meta), path(output_map_name), path(output_ped_name), emit: validated_plink
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/validate_plink.py \\
        --input-map ${map} \\
        --input-ped ${ped} \\
        --output-map ${output_map_name} \\
        --output-ped ${output_ped_name} \\
        --report ${report_name} \\
        ${args}
    """
}
