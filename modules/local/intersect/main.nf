nextflow.enable.dsl = 2

process INTERSECT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_bedtools ?: 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'}"

    input:
    tuple val(meta), path(a), path(b)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: raw

    script:
    """
    bedtools intersect -wa -wb -a ${a} -b ${b} > ${output_name}
    """
}
