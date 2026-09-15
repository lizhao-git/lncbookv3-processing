nextflow.enable.dsl = 2

process CONVERT_ALIGNMENT_FORMAT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_samtools ?: 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'}"

    input:
    tuple val(meta), path(alignment)
    val to_format
    val sort_records
    val index_output
    val threads
    val output_name

    output:
    tuple val(meta), path(output_name), emit: alignment
    tuple val(meta), path("${output_name}.bai"), optional: true, emit: index

    script:
    def sort_arg = sort_records ? '--sort' : ''
    def index_arg = index_output ? '--index' : ''
    """
    python3 ${projectDir}/scripts/format_convert/convert_alignment_format.py \\
        --input-alignment ${alignment} \\
        --output-alignment ${output_name} \\
        --to ${to_format} \\
        --threads ${threads} \\
        ${sort_arg} ${index_arg}
    """
}
