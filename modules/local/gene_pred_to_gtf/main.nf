nextflow.enable.dsl = 2

process GENE_PRED_TO_GTF {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val output_name
    val source

    output:
    tuple val(meta), path(output_name), emit: gtf

    script:
    def args = task.ext.args ?: ''
    def src = source ?: 'file'
    """
    genePredToGtf ${args} ${src} ${input_file} ${output_name}
    """
}
