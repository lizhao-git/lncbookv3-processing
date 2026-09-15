nextflow.enable.dsl = 2

process MAKE_UCSC_TRACK {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(bed)
    val track_name
    val description
    val track_type
    val color
    val visibility
    val priority
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: bed
    tuple val(meta), path(report_name), emit: report

    script:
    def type_arg = track_type ? "--track-type '${track_type}'" : ''
    def color_arg = color ? "--color '${color}'" : ''
    def visibility_arg = visibility ? "--visibility '${visibility}'" : ''
    def priority_arg = priority == null ? '' : "--priority ${priority}"
    def report_arg = report_name ? "--report ${report_name}" : ''
    """
    python3 ${projectDir}/scripts/format_convert/make_ucsc_track.py \\
        --input-bed ${bed} \\
        --output-bed ${output_name} \\
        --track-name '${track_name}' \\
        --description '${description}' \\
        ${type_arg} ${color_arg} ${visibility_arg} ${priority_arg} ${report_arg}
    """
}
