nextflow.enable.dsl = 2

process PSL2BED {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_bedops ?: 'quay.io/biocontainers/bedops:2.4.42--hd6d6fdc_1'}"

    input:
    tuple val(meta), path(input_file)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: bed

    script:
    def args = task.ext.args ?: ''
    """
    psl2bed ${args} < ${input_file} > ${output_name}
    """
}
