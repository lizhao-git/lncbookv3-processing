nextflow.enable.dsl = 2

// SMPROT: annotate SmProt functional peptide/element records onto annotation
// feature intervals.
//
// Inputs
//   records_ch    : tuple val(meta), path(smprot.tsv[.gz]) (one or more files)
//   annotation_ch : tuple val(meta), path(annotation.gtf|.gff3|.bed[.gz])
//
// The annotation is accepted in GTF, GFF3 or BED format (auto-detected from
// the file extension by ANNOTATION_FEATURES) and converted to a validated
// 10-column feature-interval BED. The SmProt side validates the TSV, filters
// the records and validates the resulting BED. Both sides are combined with
// bedtools intersect and the hits are formatted into a per-record x
// per-feature annotation table with a transcript mapping report, summarized
// by feature type, and exported as SQL / UCSC track files. Step reports are
// aggregated into pipeline_run.log and the Markdown/HTML debug report via the
// shared PIPELINE_REPORTS subworkflow.

include { ANNOTATION_FEATURES } from '../annotation_features/main'
include { PIPELINE_REPORTS } from '../pipeline_reports/main'
include { VALIDATE_SMPROT_TSV } from '../../../modules/local/validate_smprot_tsv/main'
include { FILTER_SMPROT_RECORDS } from '../../../modules/local/filter_smprot_records/main'
include { VALIDATE_BED } from '../../../modules/local/validate_bed/main'
include { INTERSECT } from '../../../modules/local/intersect/main'
include { FORMAT_SMPROT_INTERSECTIONS } from '../../../modules/local/format_smprot_intersections/main'
include { SUMMARIZE_SMPROT } from '../../../modules/local/summarize_smprot/main'
include { EXPORT_SMPROT_FORMATS } from '../../../modules/local/export_smprot_formats/main'

workflow SMPROT {
    take:
    records_ch
    annotation_ch

    main:
    def feature_types = ['gene', 'transcript', 'exon', 'intron', 'cds', 'utr']

    // --- Annotation side: validate + convert (gtf/gff3/bed) to feature intervals ---
    ANNOTATION_FEATURES(annotation_ch, feature_types)

    // --- Record side: validate the SmProt TSV, filter, validate the BED ---
    VALIDATE_SMPROT_TSV(records_ch, 'validated_smprot.tsv', 'smprot_validation_report.tsv')
    FILTER_SMPROT_RECORDS(VALIDATE_SMPROT_TSV.out.validated)
    VALIDATE_BED(FILTER_SMPROT_RECORDS.out.bed, 4, 'validated_smprot_filtered.bed', 'smprot_bed_validation_report.tsv')

    // --- Interval operation: every feature whose span overlaps a record ---
    variant_bed_out = VALIDATE_BED.out.validated_bed
    features_out = ANNOTATION_FEATURES.out.features
    intersect_in = variant_bed_out
        .combine(features_out)
        .map { meta, bed, fmeta, features -> tuple(meta, features, bed) }
    INTERSECT(intersect_in, 'smprot_intersections.tsv')

    // --- Format conversion, summary and export ---
    FORMAT_SMPROT_INTERSECTIONS(INTERSECT.out.raw, 'smprot')
    SUMMARIZE_SMPROT(FORMAT_SMPROT_INTERSECTIONS.out.annotations, 'smprot_summary.tsv')
    EXPORT_SMPROT_FORMATS(FORMAT_SMPROT_INTERSECTIONS.out.annotations, 'smprot')

    // --- Run log + debug report via the shared reporting subworkflow ---
    // NOTE: assign mixed channels to local variables before further chaining.
    reports_mix = ANNOTATION_FEATURES.out.validation_reports
        .mix(
            VALIDATE_SMPROT_TSV.out.report,
            VALIDATE_BED.out.report,
            FILTER_SMPROT_RECORDS.out.report,
            FORMAT_SMPROT_INTERSECTIONS.out.report,
            SUMMARIZE_SMPROT.out.summary
        )
    PIPELINE_REPORTS(reports_mix, FORMAT_SMPROT_INTERSECTIONS.out.annotations, 'smprot', "SmProt record annotation debug report")

    emit:
    features = ANNOTATION_FEATURES.out.features
    annotation_validation_reports = ANNOTATION_FEATURES.out.validation_reports
    validated_tsv = VALIDATE_SMPROT_TSV.out.validated
    validation_report = VALIDATE_SMPROT_TSV.out.report
    filtered_bed = VALIDATE_BED.out.validated_bed
    filtered_tsv = FILTER_SMPROT_RECORDS.out.tsv
    filter_report = FILTER_SMPROT_RECORDS.out.report
    raw_intersections = INTERSECT.out.raw
    annotations = FORMAT_SMPROT_INTERSECTIONS.out.annotations
    mapping_report = FORMAT_SMPROT_INTERSECTIONS.out.report
    summary = SUMMARIZE_SMPROT.out.summary
    sql = EXPORT_SMPROT_FORMATS.out.sql
    ucsc = EXPORT_SMPROT_FORMATS.out.ucsc
    pipeline_log = PIPELINE_REPORTS.out.pipeline_log
    debug_report_md = PIPELINE_REPORTS.out.debug_report_md
    debug_report_html = PIPELINE_REPORTS.out.debug_report_html
}
