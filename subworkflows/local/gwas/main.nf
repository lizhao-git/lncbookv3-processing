nextflow.enable.dsl = 2
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
    features_ch

    main:

    VALIDATE_GWAS_TSV(records_ch, 'validated_gwas.tsv', 'gwas_validation_report.tsv')
    FILTER_GWAS_SIGNIFICANT(VALIDATE_GWAS_TSV.out.validated)
    VALIDATE_BED(FILTER_GWAS_SIGNIFICANT.out.bed, 4, 'gwas_catalog_significant.bed', 'gwas_bed_validation_report.tsv')
    intersect_in = VALIDATE_BED.out.validated_bed.combine(features_ch).map { meta, bed, fmeta, features -> tuple(meta, features, bed) }
    INTERSECT(intersect_in, 'gwas_intersections.tsv')
    FORMAT_GWAS_INTERSECTIONS(INTERSECT.out.raw)
    SUMMARIZE_GWAS(FORMAT_GWAS_INTERSECTIONS.out.annotations, 'gtf_gwas_summary.tsv')
    EXPORT_GWAS_FORMATS(FORMAT_GWAS_INTERSECTIONS.out.annotations)


    emit:
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
}
