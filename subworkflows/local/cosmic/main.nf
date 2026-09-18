nextflow.enable.dsl = 2

// COSMIC: annotate COSMIC coding-cancer mutations onto annotation feature intervals.
//
// Inputs
//   records_ch    : tuple val(meta), path(cosmic.tsv[.gz])
//   annotation_ch : tuple val(meta), path(annotation.gtf|.gff3|.bed[.gz])
//
// The annotation is accepted in GTF, GFF3 or BED format (auto-detected from
// the file extension by ANNOTATION_FEATURES) and converted to a validated
// 10-column feature-interval BED. The COSMIC side validates the TSV, filters
// to records with pathogenic FATHMM-MKL predictions and validates the
// resulting variant BED. Both sides are combined with bedtools intersect and
// the hits are formatted into a per-mutation x per-feature annotation table,
// summarized by feature type, and exported as SQL / UCSC track files. Every
// step writes a metric/value TSV report; AGGREGATE_PIPELINE_LOGS concatenates
// them into a single pipeline_run.log and DEBUG_REPORT renders a Markdown +
// HTML debug overview.

include { ANNOTATION_FEATURES } from '../annotation_features/main'
include { PIPELINE_REPORTS } from '../pipeline_reports/main'
include { VALIDATE_COSMIC_TSV } from '../../../modules/local/validate_cosmic_tsv/main'
include { FILTER_COSMIC_FATHMM } from '../../../modules/local/filter_cosmic_fathmm/main'
include { VALIDATE_BED } from '../../../modules/local/validate_bed/main'
include { INTERSECT } from '../../../modules/local/intersect/main'
include { FORMAT_COSMIC_INTERSECTIONS } from '../../../modules/local/format_cosmic_intersections/main'
include { SUMMARIZE_ANNOTATIONS } from '../../../modules/local/summarize_annotations/main'
include { EXPORT_COSMIC_FORMATS } from '../../../modules/local/export_cosmic_formats/main'

workflow COSMIC {
    take:
    records_ch
    annotation_ch

    main:
    def feature_types = ['gene', 'transcript', 'exon', 'intron', 'cds', 'utr']

    // --- Annotation side: validate + convert (gtf/gff3/bed) to feature intervals ---
    ANNOTATION_FEATURES(annotation_ch, feature_types)

    // --- Variant side: validate the COSMIC TSV, keep pathogenic FATHMM calls ---
    VALIDATE_COSMIC_TSV(records_ch, 'validated_cosmic.tsv', 'cosmic_validation_report.tsv')
    FILTER_COSMIC_FATHMM(VALIDATE_COSMIC_TSV.out.validated)
    VALIDATE_BED(FILTER_COSMIC_FATHMM.out.bed, 4, 'validated_cosmic_pathogenic.bed', 'cosmic_bed_validation_report.tsv')

    // --- Interval operation: every feature whose span overlaps a mutation ---
    variant_bed_out = VALIDATE_BED.out.validated_bed
    features_out = ANNOTATION_FEATURES.out.features
    intersect_in = variant_bed_out
        .combine(features_out)
        .map { meta, bed, fmeta, features -> tuple(meta, features, bed) }
    INTERSECT(intersect_in, 'cosmic_intersections.tsv')

    // --- Format conversion, summary and export ---
    FORMAT_COSMIC_INTERSECTIONS(INTERSECT.out.raw, 'cosmic')
    SUMMARIZE_ANNOTATIONS(FORMAT_COSMIC_INTERSECTIONS.out.annotations, 'cosmic_summary.tsv')
    EXPORT_COSMIC_FORMATS(FORMAT_COSMIC_INTERSECTIONS.out.annotations, 'cosmic')

    // --- Run log + debug report via the shared reporting subworkflow ---
    // NOTE: assign mixed channels to local variables before further chaining.
    reports_mix = ANNOTATION_FEATURES.out.validation_reports
        .mix(
            VALIDATE_COSMIC_TSV.out.report,
            VALIDATE_BED.out.report,
            FILTER_COSMIC_FATHMM.out.report,
            SUMMARIZE_ANNOTATIONS.out.summary
        )
    PIPELINE_REPORTS(reports_mix, FORMAT_COSMIC_INTERSECTIONS.out.annotations, 'cosmic', "COSMIC mutation annotation debug report")

    emit:
    features = ANNOTATION_FEATURES.out.features
    annotation_validation_reports = ANNOTATION_FEATURES.out.validation_reports
    validated_tsv = VALIDATE_COSMIC_TSV.out.validated
    validation_report = VALIDATE_COSMIC_TSV.out.report
    filtered_bed = VALIDATE_BED.out.validated_bed
    filtered_tsv = FILTER_COSMIC_FATHMM.out.tsv
    filter_report = FILTER_COSMIC_FATHMM.out.report
    raw_intersections = INTERSECT.out.raw
    annotations = FORMAT_COSMIC_INTERSECTIONS.out.annotations
    summary = SUMMARIZE_ANNOTATIONS.out.summary
    sql = EXPORT_COSMIC_FORMATS.out.sql
    ucsc = EXPORT_COSMIC_FORMATS.out.ucsc
    pipeline_log = PIPELINE_REPORTS.out.pipeline_log
    debug_report_md = PIPELINE_REPORTS.out.debug_report_md
    debug_report_html = PIPELINE_REPORTS.out.debug_report_html
}
