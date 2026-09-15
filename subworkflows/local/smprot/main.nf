nextflow.enable.dsl = 2
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
    features_ch

    main:

    VALIDATE_SMPROT_TSV(records_ch, 'validated_smprot.tsv', 'smprot_validation_report.tsv')
    FILTER_SMPROT_RECORDS(VALIDATE_SMPROT_TSV.out.validated)
    VALIDATE_BED(FILTER_SMPROT_RECORDS.out.bed, 4, 'smprot_filtered.bed', 'smprot_bed_validation_report.tsv')
    intersect_in = VALIDATE_BED.out.validated_bed.combine(features_ch).map { meta, bed, fmeta, features -> tuple(meta, features, bed) }
    INTERSECT(intersect_in, 'smprot_intersections.tsv')
    FORMAT_SMPROT_INTERSECTIONS(INTERSECT.out.raw)
    SUMMARIZE_SMPROT(FORMAT_SMPROT_INTERSECTIONS.out.annotations, 'gtf_smprot_summary.tsv')
    EXPORT_SMPROT_FORMATS(FORMAT_SMPROT_INTERSECTIONS.out.annotations)


    emit:
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
}
