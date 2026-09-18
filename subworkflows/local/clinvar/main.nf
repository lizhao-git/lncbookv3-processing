nextflow.enable.dsl = 2

// CLINVAR: annotate ClinVar variants onto annotation feature intervals.
//
// Inputs
//   clinvar_vcf_ch : tuple val(meta), path(clinvar.vcf[.gz])
//   annotation_ch  : tuple val(meta), path(annotation.gtf|.gff3|.bed[.gz])
//
// The annotation is accepted in GTF, GFF3 or BED format (auto-detected from
// the file extension by ANNOTATION_FEATURES) and converted to a validated
// 10-column feature-interval BED. The ClinVar side validates the VCF, filters
// the records down to definite clinical-significance labels (Pathogenic,
// Benign, ...) and validates the resulting variant BED. Both sides are then
// combined with bedtools intersect and the hits are formatted into a
// per-variant x per-feature annotation table, summarized by feature type and
// clinical label, and exported as SQL / UCSC track files. Every step writes a
// metric/value TSV report; AGGREGATE_PIPELINE_LOGS concatenates them into a
// single pipeline_run.log.

include { ANNOTATION_FEATURES } from '../annotation_features/main'
include { PIPELINE_REPORTS } from '../pipeline_reports/main'
include { VALIDATE_CLINVAR_VCF } from '../../../modules/local/validate_clinvar_vcf/main'
include { FILTER_CLINVAR_LABELS } from '../../../modules/local/filter_clinvar_labels/main'
include { VALIDATE_BED } from '../../../modules/local/validate_bed/main'
include { INTERSECT } from '../../../modules/local/intersect/main'
include { FORMAT_CLINVAR_INTERSECTIONS } from '../../../modules/local/format_clinvar_intersections/main'
include { SUMMARIZE_ANNOTATIONS } from '../../../modules/local/summarize_annotations/main'
include { EXPORT_CLINVAR_FORMATS } from '../../../modules/local/export_clinvar_formats/main'

workflow CLINVAR {
    take:
    clinvar_vcf_ch
    annotation_ch

    main:
    def feature_types = ['gene', 'transcript', 'exon', 'intron', 'cds', 'utr']

    // --- Annotation side: validate + convert (gtf/gff3/bed) to feature intervals ---
    ANNOTATION_FEATURES(annotation_ch, feature_types)

    // --- Variant side: validate the ClinVar VCF, then split by clinical significance ---
    VALIDATE_CLINVAR_VCF(clinvar_vcf_ch, 'validated_clinvar.vcf', 'clinvar_validation_report.tsv')
    FILTER_CLINVAR_LABELS(VALIDATE_CLINVAR_VCF.out.validated)
    VALIDATE_BED(FILTER_CLINVAR_LABELS.out.bed, 4, 'validated_clinvar_definite_labels.bed', 'clinvar_bed_validation_report.tsv')

    // --- Interval operation: every feature whose span overlaps a variant ---
    variant_bed_out = VALIDATE_BED.out.validated_bed
    features_out = ANNOTATION_FEATURES.out.features
    intersect_in = variant_bed_out
        .combine(features_out)
        .map { meta, bed, fmeta, features -> tuple(meta, features, bed) }
    INTERSECT(intersect_in, 'clinvar_intersections.tsv')

    // --- Format conversion, summary and export ---
    FORMAT_CLINVAR_INTERSECTIONS(INTERSECT.out.raw, 'clinvar')
    SUMMARIZE_ANNOTATIONS(FORMAT_CLINVAR_INTERSECTIONS.out.annotations, 'clinvar_summary.tsv')
    EXPORT_CLINVAR_FORMATS(FORMAT_CLINVAR_INTERSECTIONS.out.annotations, 'clinvar')

    // --- Run log + debug report via the shared reporting subworkflow ---
    // Report tuples carry different meta ids (annotation basename vs VCF
    // basename); PIPELINE_REPORTS drops the meta before aggregating. The
    // annotation table is included so the debug report doubles as a data
    // spot-check.
    // NOTE: assign mixed channels to local variables before further chaining.
    reports_mix = ANNOTATION_FEATURES.out.validation_reports
        .mix(
            VALIDATE_CLINVAR_VCF.out.report,
            VALIDATE_BED.out.report,
            FILTER_CLINVAR_LABELS.out.report,
            SUMMARIZE_ANNOTATIONS.out.summary
        )
    PIPELINE_REPORTS(reports_mix, FORMAT_CLINVAR_INTERSECTIONS.out.annotations, 'clinvar', "ClinVar variant annotation debug report")

    emit:
    features = ANNOTATION_FEATURES.out.features
    annotation_validation_reports = ANNOTATION_FEATURES.out.validation_reports
    validated_vcf = VALIDATE_CLINVAR_VCF.out.validated
    validation_report = VALIDATE_CLINVAR_VCF.out.report
    filtered_bed = VALIDATE_BED.out.validated_bed
    filtered_vcf = FILTER_CLINVAR_LABELS.out.vcf
    filter_report = FILTER_CLINVAR_LABELS.out.report
    raw_intersections = INTERSECT.out.raw
    annotations = FORMAT_CLINVAR_INTERSECTIONS.out.annotations
    summary = SUMMARIZE_ANNOTATIONS.out.summary
    sql = EXPORT_CLINVAR_FORMATS.out.sql
    ucsc = EXPORT_CLINVAR_FORMATS.out.ucsc
    pipeline_log = PIPELINE_REPORTS.out.pipeline_log
    debug_report_md = PIPELINE_REPORTS.out.debug_report_md
    debug_report_html = PIPELINE_REPORTS.out.debug_report_html
}
