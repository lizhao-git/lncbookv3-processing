nextflow.enable.dsl = 2

process VALIDATE_GWAS_TSV {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(input_file)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated
    tuple val(meta), path(report_name), emit: report

    script:
    """
    python3 ${projectDir}/scripts/variant_pipeline/validate_gwas_catalog_tsv.py \
        --input-tsv ${input_file} \
        --output-tsv ${output_name} \
        --report ${report_name}
    """
}
