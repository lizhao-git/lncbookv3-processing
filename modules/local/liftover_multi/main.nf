nextflow.enable.dsl = 2

process LIFTOVER_MULTI {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${params.container_kent ?: 'lncbookv3-kent:latest'}"

    input:
    tuple val(meta), path(manifest)
    path data_root

    output:
    tuple val(meta), path('liftover_outputs'), emit: outputs_dir
    tuple val(meta), path('liftover_report.tsv'), emit: report

    script:
    def root_arg = data_root ? "--data-root ${data_root}" : ''
    """
    python3 ${projectDir}/scripts/format_convert/liftover_multi.py \\
        --manifest ${manifest} \\
        ${root_arg} \\
        --output-dir liftover_outputs \\
        --report liftover_report.tsv
    """
}
