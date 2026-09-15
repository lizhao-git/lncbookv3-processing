nextflow.enable.dsl = 2
include { VALIDATE_CLINVAR_VCF } from '../../../modules/local/validate_clinvar_vcf/main'
include { FILTER_CLINVAR_LABELS } from '../../../modules/local/filter_clinvar_labels/main'
include { VALIDATE_BED } from '../../../modules/local/validate_bed/main'
include { INTERSECT } from '../../../modules/local/intersect/main'
include { FORMAT_CLINVAR_INTERSECTIONS } from '../../../modules/local/format_clinvar_intersections/main'
include { SUMMARIZE_ANNOTATIONS } from '../../../modules/local/summarize_annotations/main'
include { EXPORT_CLINVAR_FORMATS } from '../../../modules/local/export_clinvar_formats/main'

workflow CLINVAR {
    take:
    records_ch
    features_ch

    main:

    VALIDATE_CLINVAR_VCF(records_ch, 'validated_clinvar.vcf', 'clinvar_validation_report.tsv')
    FILTER_CLINVAR_LABELS(VALIDATE_CLINVAR_VCF.out.validated)
    VALIDATE_BED(FILTER_CLINVAR_LABELS.out.bed, 4, 'clinvar_definite_labels.bed', 'clinvar_bed_validation_report.tsv')
    intersect_in = VALIDATE_BED.out.validated_bed.combine(features_ch).map { meta, bed, fmeta, features -> tuple(meta, features, bed) }
    INTERSECT(intersect_in, 'clinvar_intersections.tsv')
    FORMAT_CLINVAR_INTERSECTIONS(INTERSECT.out.raw)
    SUMMARIZE_ANNOTATIONS(FORMAT_CLINVAR_INTERSECTIONS.out.annotations, 'gtf_clinvar_summary.tsv')
    EXPORT_CLINVAR_FORMATS(FORMAT_CLINVAR_INTERSECTIONS.out.annotations)


    emit:
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
}
