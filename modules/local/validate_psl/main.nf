nextflow.enable.dsl = 2

process VALIDATE_PSL {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(psl), path(target_sizes)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated_psl
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    def sizes_arg = target_sizes ? "--target-sizes ${target_sizes}" : ''
    """
    python3 ${projectDir}/scripts/format_convert/validate_psl.py \\
        --input-psl ${psl} \\
        --output-psl ${output_name} \\
        --report ${report_name} \\
        ${sizes_arg} \\
        ${args}
    """
}
