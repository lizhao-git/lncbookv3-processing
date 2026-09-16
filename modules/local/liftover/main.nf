nextflow.enable.dsl = 2

process LIFTOVER {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_track)
    path chain
    val min_match
    val output_name
    val unmapped_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: lifted_bed
    tuple val(meta), path(unmapped_name), emit: unmapped_bed
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/liftover.py \\
        --input-track ${input_track} \\
        --chain ${chain} \\
        --output-bed ${output_name} \\
        --unmapped-bed ${unmapped_name} \\
        --report ${report_name} \\
        --min-match ${min_match} \\
        ${args}
    """
}
