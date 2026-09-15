nextflow.enable.dsl = 2

process VCF_SPLIT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_python ?: 'python:3.12-slim'}"

    input:
    tuple val(meta), path(vcf)
    val mode
    val prefix
    val out_dir
    val report_name

    output:
    tuple val(meta), path(out_dir), emit: split_files
    tuple val(meta), path(report_name), emit: report

    script:
    """
    python3 ${projectDir}/scripts/format_convert/vcf_split.py \\
        --input-vcf ${vcf} \\
        --out-dir ${out_dir} \\
        --mode ${mode} \\
        --prefix ${prefix} \\
        --report ${report_name}
    """
}
