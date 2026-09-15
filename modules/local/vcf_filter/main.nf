nextflow.enable.dsl = 2

process VCF_FILTER {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(vcf)
    val regions
    val min_qual
    val min_dp
    val keep_types
    val require_pass
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: vcf
    tuple val(meta), path(report_name), emit: report

    script:
    def region_args = regions ? regions.collect { "--region '${it}'" }.join(' ') : ''
    def type_args = keep_types ? keep_types.collect { "--keep-type ${it}" }.join(' ') : ''
    def qual_arg = min_qual == null ? '' : "--min-qual ${min_qual}"
    def dp_arg = min_dp == null ? '' : "--min-dp ${min_dp}"
    def pass_arg = require_pass ? '--require-pass' : ''
    """
    python3 ${projectDir}/scripts/format_convert/vcf_filter.py \\
        --input-vcf ${vcf} \\
        --output-vcf ${output_name} \\
        --report ${report_name} \\
        ${region_args} ${qual_arg} ${dp_arg} ${type_args} ${pass_arg}
    """
}
