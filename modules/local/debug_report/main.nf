nextflow.enable.dsl = 2

process DEBUG_REPORT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(reports, stageAs: 'inputs/*')
    val md_name
    val html_name
    val report_title

    output:
    tuple val(meta), path(md_name), emit: report_md
    tuple val(meta), path(html_name), emit: report_html

    script:
    def report_args = reports.collect { " '--input' '${it}'" }.join(' \\\n        ')
    """
    python3 ${projectDir}/scripts/variant_pipeline/generate_debug_report.py \\
        --title '${report_title}' \\
        --output-md ${md_name} \\
        --output-html ${html_name} \\
        ${report_args}
    """
}
