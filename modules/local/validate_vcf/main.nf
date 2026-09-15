nextflow.enable.dsl = 2

process VALIDATE_VCF {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(vcf)
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: validated_vcf
    tuple val(meta), path(report_name), emit: report

    script:
    """
    python3 ${projectDir}/scripts/format_convert/validate_vcf.py \\
        --input-vcf ${vcf} \\
        --output-vcf ${output_name} \\
        --report ${report_name}
    """
}
