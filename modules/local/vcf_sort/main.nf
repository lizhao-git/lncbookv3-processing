nextflow.enable.dsl = 2

process VCF_SORT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(vcf)
    val contig_order
    path fasta_index
    val output_name
    val report_name

    output:
    tuple val(meta), path(output_name), emit: vcf
    tuple val(meta), path(report_name), emit: report

    script:
    def order_arg = contig_order ? "--contig-order '${contig_order}'" : ''
    def fai_arg = fasta_index ? "--fasta-index ${fasta_index}" : ''
    def report_arg = report_name ? "--report ${report_name}" : ''
    """
    python3 ${projectDir}/scripts/format_convert/vcf_sort.py \\
        --input-vcf ${vcf} \\
        --output-vcf ${output_name} \\
        ${order_arg} ${fai_arg} ${report_arg}
    """
}
