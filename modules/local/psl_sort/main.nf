nextflow.enable.dsl = 2

process PSL_SORT {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: psl

    script:
    def args = task.ext.args ?: ''
    def dirs = task.ext.dirs ?: 'dirs'
    """
    mkdir -p tmp_pslsort
    pslSort ${dirs} ${args} ${output_name} tmp_pslsort ${input_file}
    rm -rf tmp_pslsort
    """
}
