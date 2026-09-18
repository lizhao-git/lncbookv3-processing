nextflow.enable.dsl = 2

process BED_TO_BEDGRAPH {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(input_file)
    path sizes
    val output_name

    output:
    tuple val(meta), path(output_name), emit: bedgraph

    script:
    def args = task.ext.args ?: ''
    def scale = task.ext.scale ?: ''
    def scale_arg = scale ? "-scale ${scale}" : ''
    """
    bedtools genomecov \\
        -i ${input_file} \\
        -g ${sizes} \\
        -bg \\
        ${scale_arg} \\
        ${args} \\
        | LC_ALL=C sort -k1,1 -k2,2n > ${output_name}
    """
}
