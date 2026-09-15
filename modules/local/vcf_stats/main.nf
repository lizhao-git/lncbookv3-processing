nextflow.enable.dsl = 2

process VCF_STATS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(vcf)
    val report_name

    output:
    tuple val(meta), path(report_name), emit: report

    script:
    """
    python3 ${projectDir}/scripts/format_convert/vcf_stats.py --input-vcf ${vcf} --report ${report_name}
    """
}
