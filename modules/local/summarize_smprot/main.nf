nextflow.enable.dsl = 2

process SUMMARIZE_SMPROT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(annotations)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: summary

    script:
    """
    python3 ${projectDir}/scripts/smprot_pipeline/summarize_smprot_annotations.py --input-tsv ${annotations} --output-summary ${output_name}
    """
}
