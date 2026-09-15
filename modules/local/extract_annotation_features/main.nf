nextflow.enable.dsl = 2

process EXTRACT_ANNOTATION_FEATURES {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(annotation)
    val format
    val features
    val output_name

    output:
    tuple val(meta), path(output_name), emit: features_bed

    script:
    def feature_args = features ? features.collect { "--feature ${it}" }.join(' ') : ''
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/scripts/format_convert/extract_annotation_features.py \\
        --input-annotation ${annotation} \\
        --output-bed ${output_name} \\
        --format ${format} \\
        ${feature_args} \\
        ${args}
    """
}
