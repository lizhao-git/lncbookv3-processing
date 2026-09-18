nextflow.enable.dsl = 2

process AGGREGATE_PIPELINE_LOGS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(reports, stageAs: 'inputs/*')
    val log_name

    output:
    tuple val(meta), path(log_name), emit: log

    script:
    // `reports` elements render as their staged paths (inputs/<name>) because
    // of the stageAs directive -- do not rebuild paths with it.name.
    def report_args = reports.collect { " '--input' '${it}'" }.join(' \\\n        ')
    def title = task.ext.title ?: 'lncbookv3-processing pipeline run log'
    """
    python3 ${projectDir}/scripts/variant_pipeline/aggregate_run_log.py \\
        --output-log ${log_name} \\
        --title '${title}' \\
        ${report_args}
    """
}
