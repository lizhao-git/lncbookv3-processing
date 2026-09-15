nextflow.enable.dsl = 2

process CERNA {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(manifest)
    path data_root

    output:
    tuple val(meta), path('cerna_outputs'), emit: outputs_dir


    script:
    def root_arg = data_root ? "--data-root ${data_root}" : ''
    """
    python3 ${projectDir}/scripts/cerna_pipeline/pipeline.py \
        --manifest ${manifest} \
        ${root_arg} \
        --output-dir cerna_outputs \

    """
}
