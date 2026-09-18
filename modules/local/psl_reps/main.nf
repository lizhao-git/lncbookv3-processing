nextflow.enable.dsl = 2

process PSL_REPS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: psl
    tuple val(meta), path(report_name), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    pslReps ${args} ${input_file} ${output_name} ${report_name}
    """
}
