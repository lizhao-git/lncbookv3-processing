nextflow.enable.dsl = 2

process BED_TO_PSL {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val output_name
    path target_sizes

    output:
    tuple val(meta), path(output_name), emit: psl

    script:
    def args = task.ext.args ?: ''
    def sizes_arg = target_sizes && target_sizes.size() > 0 ? "--target-sizes ${target_sizes}" : ""
    """
    python3 ${projectDir}/scripts/format_convert/bed_to_psl.py \\
        --input-bed ${input_file} \\
        --output-psl ${output_name} \\
        ${sizes_arg} \\
        ${args}
    """
}
