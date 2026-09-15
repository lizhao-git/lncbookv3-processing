nextflow.enable.dsl = 2

process BEDGRAPH_TO_BIGWIG {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(input_file)
    path chrom_sizes
    val clip
    val report_name
    val output_name

    output:
    tuple val(meta), path(output_name), emit: output
    tuple val(meta), path(report_name), optional: true, emit: report

    script:
    def clip_arg = clip ? '--clip' : ''; def report_arg = report_name ? "--report ${report_name}" : ''
    """
    python3 ${projectDir}/scripts/format_convert/bedgraph_to_bigwig.py \
        --input-bedgraph ${input_file} \
        --output-bigwig ${output_name} \
        --chrom-sizes ${chrom_sizes} \
        ${clip_arg} ${report_arg}
    """
}
