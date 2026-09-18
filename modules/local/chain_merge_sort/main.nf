nextflow.enable.dsl = 2

process CHAIN_MERGE_SORT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: chain

    script:
    def args = task.ext.args ?: ''
    """
    chainMergeSort ${args} ${input_file} > ${output_name}
    """
}
