nextflow.enable.dsl = 2

process BED_TO_BED12 {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: bed12

    script:
    def args = task.ext.args ?: ''
    def rgb = task.ext.rgb ?: '0'
    """
    python3 ${projectDir}/scripts/format_convert/bed_to_bed12.py \\
        --input-bed ${input_file} \\
        --output-bed ${output_name} \\
        --rgb ${rgb} \\
        ${args}
    """
}
