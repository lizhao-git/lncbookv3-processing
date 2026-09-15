nextflow.enable.dsl = 2

process BIGTRACK_TO_TEXT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    val report_name
    val output_name

    output:
    tuple val(meta), path(output_name), emit: output
    tuple val(meta), path(report_name), optional: true, emit: report

    script:
    def report_arg = report_name ? "--report ${report_name}" : ''
    """
    python3 ${projectDir}/scripts/format_convert/bigtrack_to_text.py \
        --input-track ${input_file} \
        --output-bed ${output_name} \
        ${report_arg}
    """
}
