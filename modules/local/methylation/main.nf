nextflow.enable.dsl = 2

process METHYLATION {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(manifest)
    path data_root

    output:
    tuple val(meta), path('methylation_outputs'), emit: outputs_dir
    tuple val(meta), path('methylation_report.tsv'), emit: report

    script:
    def root_arg = data_root ? "--data-root ${data_root}" : ''
    """
    python3 ${projectDir}/scripts/methylation/pipeline.py \
        --manifest ${manifest} \
        ${root_arg} \
        --output-dir methylation_outputs \
        --report methylation_report.tsv
    """
}
