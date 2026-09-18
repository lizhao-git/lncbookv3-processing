nextflow.enable.dsl = 2

// PIPELINE_REPORTS: shared run-log + debug-report composition.
//
// Extracts the reporting tail that every analysis branch (clinvar, cosmic,
// gwas, smprot, methylation, conservation) shares:
//
//   1. AGGREGATE_PIPELINE_LOGS concatenates every step report (metric/value
//      TSVs and summaries) into a single `pipeline_run.log`, one
//      `===== <report name> =====` section per step.
//   2. DEBUG_REPORT renders the same reports plus optional data tables into
//      `pipeline_debug_report.md` / `.html` — overview table with OK/WARN/FAIL
//      badges per step, errors and warnings listed first, per-step details.
//
// Inputs
//   reports_in : channel of tuple val(meta), path(report.tsv) — every step
//                report / summary to include. Meta ids may differ between
//                branches (annotation basename vs records basename); they are
//                dropped before aggregation.
//   data_in    : channel of tuple val(meta), path(data.tsv) — optional data
//                tables (e.g. the annotation TSV) included in the debug report
//                for spot-checking. Pass Channel.empty() when there is none.
//   branch_id  : val — meta id used for the aggregated tuple and the publish
//                subdirectory.
//   title      : val — title line of the debug report.
//
// Emits: pipeline_log, debug_report_md, debug_report_html (each a single
// tuple val(meta), path(file) item per invocation).

include { AGGREGATE_PIPELINE_LOGS } from '../../../modules/local/aggregate_pipeline_logs/main'
include { DEBUG_REPORT } from '../../../modules/local/debug_report/main'

workflow PIPELINE_REPORTS {
    take:
    reports_in
    data_in
    branch_id
    title

    main:
    // Flatten the mixed report tuples (different meta ids across branches)
    // into plain paths, collect them, and sort by file name so the run log
    // sections appear in a stable order.
    report_list = reports_in
        .map { meta, file -> file }
        .toList()
        .map { files -> files.sort { it.name } }

    // Optional data tables: `toList()` on an empty channel never emits, so
    // `ifEmpty([])` guarantees a value and the downstream merge always runs.
    data_list = data_in
        .map { meta, file -> file }
        .toList()
        .ifEmpty([])
        .map { files -> files.sort { it.name } }

    AGGREGATE_PIPELINE_LOGS(report_list.map { files -> tuple([id: branch_id], files) }, 'pipeline_run.log')

    // Merge report + data lists into one flat, sorted list for DEBUG_REPORT.
    // NOTE: chain operators on local variables, never directly on
    // `WORKFLOW.out.*` — the DSL cannot resolve that form.
    report_list_out = report_list
    data_list_out = data_list
    all_files = report_list_out
        .mix(data_list_out)
        .toList()
        .map { lists -> lists.flatten().sort { it.name } }
        .map { files -> tuple([id: branch_id], files) }
    DEBUG_REPORT(all_files, 'pipeline_debug_report.md', 'pipeline_debug_report.html', title)

    emit:
    pipeline_log = AGGREGATE_PIPELINE_LOGS.out.log
    debug_report_md = DEBUG_REPORT.out.report_md
    debug_report_html = DEBUG_REPORT.out.report_html
}
