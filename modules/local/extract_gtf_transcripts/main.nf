nextflow.enable.dsl = 2

process EXTRACT_GTF_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(gtf)
    val output_name

    output:
    tuple val(meta), path(output_name), emit: transcripts_bed

    script:
    """
    python3 ${projectDir}/scripts/format_convert/extract_gtf_transcripts.py \\
        --input-gtf ${gtf} \\
        --output-bed ${output_name}
    """
}
