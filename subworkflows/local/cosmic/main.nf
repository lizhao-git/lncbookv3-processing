nextflow.enable.dsl = 2
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
    features_ch

    main:

    VALIDATE_COSMIC_TSV(records_ch, 'validated_cosmic.tsv', 'cosmic_validation_report.tsv')
    FILTER_COSMIC_FATHMM(VALIDATE_COSMIC_TSV.out.validated)
    VALIDATE_BED(FILTER_COSMIC_FATHMM.out.bed, 4, 'cosmic_pathogenic.bed', 'cosmic_bed_validation_report.tsv')
    intersect_in = VALIDATE_BED.out.validated_bed.combine(features_ch).map { meta, bed, fmeta, features -> tuple(meta, features, bed) }
    INTERSECT(intersect_in, 'cosmic_intersections.tsv')
    FORMAT_COSMIC_INTERSECTIONS(INTERSECT.out.raw)
    SUMMARIZE_ANNOTATIONS(FORMAT_COSMIC_INTERSECTIONS.out.annotations, 'gtf_cosmic_summary.tsv')
    EXPORT_COSMIC_FORMATS(FORMAT_COSMIC_INTERSECTIONS.out.annotations)


    emit:
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
}
