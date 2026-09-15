nextflow.enable.dsl = 2

process PSLMAP_CHAIN {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(query_psl), path(chain)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: psl

    script:
    def args = task.ext.args ?: ''
    """
    pslMap ${args} ${query_psl} ${chain} ${output_name}
    """
}
