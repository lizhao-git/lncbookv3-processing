nextflow.enable.dsl = 2

process BED_TO_BIGBED {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    path chrom_sizes
    val bed_type
    val report_name
    val output_name

    output:
    tuple val(meta), path(output_name), emit: output
    tuple val(meta), path(report_name), optional: true, emit: report

    script:
    def type_arg = bed_type ? "--bed-type ${bed_type}" : ''; def report_arg = report_name ? "--report ${report_name}" : ''
    """
    python3 ${projectDir}/scripts/format_convert/bed_to_bigbed.py \
        --input-bed ${input_file} \
        --output-bigbed ${output_name} \
        --chrom-sizes ${chrom_sizes} \
        ${type_arg} ${report_arg}
    """
}
