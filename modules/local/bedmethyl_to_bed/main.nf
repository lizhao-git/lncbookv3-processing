nextflow.enable.dsl = 2

process BEDMETHYL_TO_BED {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: bed

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/bedmethyl_to_bed.py \\
        --input-bedmethyl ${input_file} \\
        --output-bed ${output_name} \\
        ${args}
    """
}
