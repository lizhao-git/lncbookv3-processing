nextflow.enable.dsl = 2

// GWAS: annotate genome-wide significant GWAS Catalog SNPs onto annotation
// feature intervals.
//
// Inputs
//   records_ch    : tuple val(meta), path(gwas.tsv[.gz])
//   annotation_ch : tuple val(meta), path(annotation.gtf|.gff3|.bed[.gz])
//
// The annotation is accepted in GTF, GFF3 or BED format (auto-detected from
// the file extension by ANNOTATION_FEATURES) and converted to a validated
// 10-column feature-interval BED. The GWAS side validates the TSV, filters to
// genome-wide significant associations and validates the resulting variant
// BED. Both sides are combined with bedtools intersect and the hits are
// formatted into a per-SNP x per-feature annotation table, summarized by
// feature type, and exported as SQL / UCSC track files. Every step writes a
// metric/value TSV report; AGGREGATE_PIPELINE_LOGS concatenates them into a
// single pipeline_run.log and DEBUG_REPORT renders a Markdown + HTML debug
// overview.

include { ANNOTATION_FEATURES } from '../annotation_features/main'
include { PIPELINE_REPORTS } from '../pipeline_reports/main'
include { VALIDATE_GWAS_TSV } from '../../../modules/local/validate_gwas_tsv/main'
include { FILTER_GWAS_SIGNIFICANT } from '../../../modules/local/filter_gwas_significant/main'
include { VALIDATE_BED } from '../../../modules/local/validate_bed/main'
include { INTERSECT } from '../../../modules/local/intersect/main'
include { FORMAT_GWAS_INTERSECTIONS } from '../../../modules/local/format_gwas_intersections/main'
include { SUMMARIZE_GWAS } from '../../../modules/local/summarize_gwas/main'
include { EXPORT_GWAS_FORMATS } from '../../../modules/local/export_gwas_formats/main'

workflow GWAS {
    take:
    records_ch
    annotation_ch

    main:
    def feature_types = ['gene', 'transcript', 'exon', 'intron', 'cds', 'utr']

    // --- Annotation side: validate + convert (gtf/gff3/bed) to feature intervals ---
    ANNOTATION_FEATURES(annotation_ch, feature_types)

    // --- Variant side: validate the GWAS Catalog TSV, keep significant hits ---
    VALIDATE_GWAS_TSV(records_ch, 'validated_gwas.tsv', 'gwas_validation_report.tsv')
    FILTER_GWAS_SIGNIFICANT(VALIDATE_GWAS_TSV.out.validated)
    VALIDATE_BED(FILTER_GWAS_SIGNIFICANT.out.bed, 4, 'validated_gwas_catalog_significant.bed', 'gwas_bed_validation_report.tsv')

    // --- Interval operation: every feature whose span overlaps a SNP ---
    variant_bed_out = VALIDATE_BED.out.validated_bed
    features_out = ANNOTATION_FEATURES.out.features
    intersect_in = variant_bed_out
        .combine(features_out)
        .map { meta, bed, fmeta, features -> tuple(meta, features, bed) }
    INTERSECT(intersect_in, 'gwas_intersections.tsv')

    // --- Format conversion, summary and export ---
    FORMAT_GWAS_INTERSECTIONS(INTERSECT.out.raw, 'gwas')
    SUMMARIZE_GWAS(FORMAT_GWAS_INTERSECTIONS.out.annotations, 'gwas_summary.tsv')
    EXPORT_GWAS_FORMATS(FORMAT_GWAS_INTERSECTIONS.out.annotations, 'gwas')

    // --- Run log + debug report via the shared reporting subworkflow ---
    // NOTE: assign mixed channels to local variables before further chaining.
    reports_mix = ANNOTATION_FEATURES.out.validation_reports
        .mix(
            VALIDATE_GWAS_TSV.out.report,
            VALIDATE_BED.out.report,
            FILTER_GWAS_SIGNIFICANT.out.report,
            SUMMARIZE_GWAS.out.summary
        )
    PIPELINE_REPORTS(reports_mix, FORMAT_GWAS_INTERSECTIONS.out.annotations, 'gwas', "GWAS Catalog SNP annotation debug report")

    emit:
    features = ANNOTATION_FEATURES.out.features
    annotation_validation_reports = ANNOTATION_FEATURES.out.validation_reports
    validated_tsv = VALIDATE_GWAS_TSV.out.validated
    validation_report = VALIDATE_GWAS_TSV.out.report
    filtered_bed = VALIDATE_BED.out.validated_bed
    filtered_tsv = FILTER_GWAS_SIGNIFICANT.out.tsv
    filter_report = FILTER_GWAS_SIGNIFICANT.out.report
    raw_intersections = INTERSECT.out.raw
    annotations = FORMAT_GWAS_INTERSECTIONS.out.annotations
    summary = SUMMARIZE_GWAS.out.summary
    sql = EXPORT_GWAS_FORMATS.out.sql
    ucsc = EXPORT_GWAS_FORMATS.out.ucsc
    pipeline_log = PIPELINE_REPORTS.out.pipeline_log
    debug_report_md = PIPELINE_REPORTS.out.debug_report_md
    debug_report_html = PIPELINE_REPORTS.out.debug_report_html
}
