nextflow.enable.dsl = 2

include { VALIDATE_ANNOTATION } from '../../../modules/local/validate_annotation/main'
include { EXTRACT_ANNOTATION_FEATURES } from '../../../modules/local/extract_annotation_features/main'
include { VALIDATE_BED } from '../../../modules/local/validate_bed/main'

workflow ANNOTATION_INTERVALS {
    take:
    annotation_ch
    format
    features

    main:
    VALIDATE_ANNOTATION(annotation_ch, format, 'validated_annotation.txt', 'annotation_validation_report.tsv')
    EXTRACT_ANNOTATION_FEATURES(VALIDATE_ANNOTATION.out.validated_annotation, format, features, 'annotation_features.bed')
    VALIDATE_BED(EXTRACT_ANNOTATION_FEATURES.out.features_bed, 10, 'validated_annotation_features.bed', 'annotation_features_bed_validation_report.tsv')

    emit:
    validated_annotation = VALIDATE_ANNOTATION.out.validated_annotation
    annotation_validation_report = VALIDATE_ANNOTATION.out.report
    features_bed = VALIDATE_BED.out.validated_bed
    features_bed_validation_report = VALIDATE_BED.out.report
}
