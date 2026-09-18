nextflow.enable.dsl = 2

process BED_TO_GTF {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: gtf

    script:
    def args = task.ext.args ?: ''
    def source = task.ext.source ?: 'bed_to_gtf'
    def feature = task.ext.feature ?: 'exon'
    """
    python3 ${projectDir}/scripts/format_convert/bed_to_gtf.py \\
        --input-bed ${input_file} \\
        --output-gtf ${output_name} \\
        --source ${source} \\
        --feature ${feature} \\
        ${args}
    """
}
