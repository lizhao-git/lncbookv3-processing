nextflow.enable.dsl = 2

process CHAIN_TO_PSL {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(chain_file), path(target_seqs), path(query_seqs)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: psl

    script:
    def args = task.ext.args ?: ''
    """
    chainToPsl ${args} ${chain_file} /dev/null /dev/null ${target_seqs} ${query_seqs} ${output_name}
    """
}
