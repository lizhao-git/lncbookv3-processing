nextflow.enable.dsl = 2

include { VALIDATE_ANNOTATION } from '../../../modules/local/validate_annotation/main'
include { EXTRACT_ANNOTATION_FEATURES } from '../../../modules/local/extract_annotation_features/main'
include { VALIDATE_BED } from '../../../modules/local/validate_bed/main'

workflow REFERENCE {
    take:
    gtf_ch

    main:
    VALIDATE_ANNOTATION(gtf_ch, 'gtf', 'validated.gtf', 'gtf_validation_report.tsv')
    EXTRACT_ANNOTATION_FEATURES(VALIDATE_ANNOTATION.out.validated_annotation, 'gtf', ['gene','transcript','exon','intron','cds','utr'], 'gtf_features.bed')
    VALIDATE_BED(EXTRACT_ANNOTATION_FEATURES.out.features_bed, 10, 'validated_gtf_features.bed', 'gtf_features_validation_report.tsv')

    emit:
    validated_gtf = VALIDATE_ANNOTATION.out.validated_annotation
    gtf_validation_report = VALIDATE_ANNOTATION.out.report
    features = VALIDATE_BED.out.validated_bed
    features_validation_report = VALIDATE_BED.out.report
}
