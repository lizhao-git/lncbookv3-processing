nextflow.enable.dsl = 2

process BED_TO_GFF3 {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(bed)
    val feature
    val source
    val output_name

    output:
    tuple val(meta), path(output_name), emit: gff3

    script:
    """
    python3 ${projectDir}/scripts/format_convert/bed_to_gff3.py \\
        --input-bed ${bed} \\
        --output-gff3 ${output_name} \\
        --feature ${feature} \\
        --source ${source}
    """
}
