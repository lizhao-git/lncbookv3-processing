nextflow.enable.dsl = 2

// CONSERVATION: JSON-manifest-driven lncRNA conservation analysis, wrapped
// with the shared run-log + debug-report tail.
//
// Inputs
//   manifest_ch  : tuple val(meta), path(manifest.json)
//   data_root_ch : path to the data root directory (or an empty value)

include { CONSERVATION as RUN_CONSERVATION } from '../../../modules/local/conservation/main'
include { PIPELINE_REPORTS } from '../pipeline_reports/main'

workflow CONSERVATION {
    take:
    manifest_ch
    data_root_ch

    main:
    RUN_CONSERVATION(manifest_ch, data_root_ch)
    PIPELINE_REPORTS(RUN_CONSERVATION.out.report, Channel.empty(), 'conservation', "LncRNA conservation analysis debug report")

    emit:
    outputs_dir = RUN_CONSERVATION.out.outputs_dir
    report = RUN_CONSERVATION.out.report
    pipeline_log = PIPELINE_REPORTS.out.pipeline_log
    debug_report_md = PIPELINE_REPORTS.out.debug_report_md
    debug_report_html = PIPELINE_REPORTS.out.debug_report_html
}
