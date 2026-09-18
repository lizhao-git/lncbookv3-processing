nextflow.enable.dsl = 2

// METHYLATION: JSON-manifest-driven methylation analysis, wrapped with the
// shared run-log + debug-report tail.
//
// Inputs
//   manifest_ch  : tuple val(meta), path(manifest.json)
//   data_root_ch : path to the data root directory (or an empty value)

include { METHYLATION as RUN_METHYLATION } from '../../../modules/local/methylation/main'
include { PIPELINE_REPORTS } from '../pipeline_reports/main'

workflow METHYLATION {
    take:
    manifest_ch
    data_root_ch

    main:
    RUN_METHYLATION(manifest_ch, data_root_ch)
    PIPELINE_REPORTS(RUN_METHYLATION.out.report, Channel.empty(), 'methylation', "Methylation analysis debug report")

    emit:
    outputs_dir = RUN_METHYLATION.out.outputs_dir
    report = RUN_METHYLATION.out.report
    pipeline_log = PIPELINE_REPORTS.out.pipeline_log
    debug_report_md = PIPELINE_REPORTS.out.debug_report_md
    debug_report_html = PIPELINE_REPORTS.out.debug_report_html
}
