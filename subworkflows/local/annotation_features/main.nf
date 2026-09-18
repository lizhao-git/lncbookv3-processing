nextflow.enable.dsl = 2

// ANNOTATION_FEATURES: accept an annotation in GTF, GFF3 or BED format and
// emit a validated 10-column feature-interval BED, regardless of input format.
//
// Inputs
//   annotation_ch : tuple val(meta), path(annotation.gtf|.gff3|.bed[.gz])
//   features      : list of feature types to extract (gtf/gff3 branch only),
//                   e.g. ['gene','transcript','exon','intron','cds','utr']
//
// Format is detected from the file extension:
//   *.gtf[.gz]                     -> gtf  branch
//   *.gff | *.gff3 [.gz]           -> gff3 branch
//   anything else (treated as BED) -> bed  branch
//
// gtf/gff3 branch : VALIDATE_ANNOTATION -> EXTRACT_ANNOTATION_FEATURES -> VALIDATE_BED(10)
// bed branch      : VALIDATE_BED(>=3)   -> BED_TO_FEATURES              -> VALIDATE_BED(10)
//
// Every step writes a metric/value TSV report; all reports are emitted so they
// can be aggregated into a pipeline run log.
//
// NOTE: each process is invoked exactly once per workflow context (nested
// helper workflows below) because referencing `PROC.out.*` from a process that
// was invoked multiple times in the same context merges all invocations'
// outputs and breaks per-branch routing.

include { VALIDATE_ANNOTATION } from '../../../modules/local/validate_annotation/main'
include { EXTRACT_ANNOTATION_FEATURES } from '../../../modules/local/extract_annotation_features/main'
include { VALIDATE_BED } from '../../../modules/local/validate_bed/main'
include { BED_TO_FEATURES } from '../../../modules/local/bed_to_features/main'

// GTF / GFF3 -> validated feature intervals
workflow FEATURES_FROM_ANNOTATION {
    take:
    ann_ch
    features

    main:
    VALIDATE_ANNOTATION(ann_ch, 'auto', 'validated_annotation.gtf', 'annotation_validation_report.tsv')
    EXTRACT_ANNOTATION_FEATURES(VALIDATE_ANNOTATION.out.validated_annotation, 'auto', features, 'extracted_features.bed')
    VALIDATE_BED(EXTRACT_ANNOTATION_FEATURES.out.features_bed, 10, 'validated_features.bed', 'features_bed_validation_report.tsv')

    emit:
    features = VALIDATE_BED.out.validated_bed
    ann_report = VALIDATE_ANNOTATION.out.report
    bed_report = VALIDATE_BED.out.report
}

// stage 1 of the BED branch: structural validation of the raw input BED
workflow BED_VALIDATE_INPUT {
    take:
    bed_ch

    main:
    VALIDATE_BED(bed_ch, 3, 'validated_input.bed', 'input_bed_validation_report.tsv')

    emit:
    bed = VALIDATE_BED.out.validated_bed
    report = VALIDATE_BED.out.report
}

// stage 2 of the BED branch: normalize to the feature schema and re-validate
workflow BED_NORMALIZE_FEATURES {
    take:
    bed_ch

    main:
    BED_TO_FEATURES(bed_ch, 'normalized_features.bed', 'bed_normalization_report.tsv')
    VALIDATE_BED(BED_TO_FEATURES.out.features_bed, 10, 'validated_features.bed', 'bed_features_validation_report.tsv')

    emit:
    features = VALIDATE_BED.out.validated_bed
    normalize_report = BED_TO_FEATURES.out.report
    bed_features_report = VALIDATE_BED.out.report
}

// plain BED -> normalized feature intervals
workflow FEATURES_FROM_BED {
    take:
    bed_ch

    main:
    BED_VALIDATE_INPUT(bed_ch)
    BED_NORMALIZE_FEATURES(BED_VALIDATE_INPUT.out.bed)

    emit:
    features = BED_NORMALIZE_FEATURES.out.features
    input_bed_report = BED_VALIDATE_INPUT.out.report
    normalize_report = BED_NORMALIZE_FEATURES.out.normalize_report
    bed_features_report = BED_NORMALIZE_FEATURES.out.bed_features_report
}

workflow ANNOTATION_FEATURES {
    take:
    annotation_ch
    features

    main:
    // --- dispatch by file extension ---
    fmt_ch = annotation_ch
        .map { meta, ann ->
            def n = ann.name.toLowerCase()
            def fmt = n.endsWith('.gtf') || n.endsWith('.gtf.gz') ? 'gtf'
                    : (n.endsWith('.gff') || n.endsWith('.gff.gz') || n.endsWith('.gff3') || n.endsWith('.gff3.gz') ? 'gff3' : 'bed')
            tuple(meta, ann, fmt)
        }

    ann_branch = fmt_ch.filter { meta, ann, fmt -> fmt != 'bed' }.map { meta, ann, fmt -> tuple(meta, ann) }
    bed_branch = fmt_ch.filter { meta, ann, fmt -> fmt == 'bed' }

    FEATURES_FROM_ANNOTATION(ann_branch, features)
    FEATURES_FROM_BED(bed_branch.map { meta, ann, fmt -> tuple(meta, ann) })

    // --- unify the two branches ---
    // NOTE: chained operators directly on `WORKFLOW.out.*` / `PROC.out.*` are
    // not resolvable by the DSL; assign each output to a local variable first.
    ann_features_out = FEATURES_FROM_ANNOTATION.out.features
    bed_features_out = FEATURES_FROM_BED.out.features
    ann_r1 = FEATURES_FROM_ANNOTATION.out.ann_report
    ann_r3 = FEATURES_FROM_ANNOTATION.out.bed_report
    bed_r1 = FEATURES_FROM_BED.out.input_bed_report
    bed_r2 = FEATURES_FROM_BED.out.normalize_report
    bed_r3 = FEATURES_FROM_BED.out.bed_features_report
    features_ch = ann_features_out.mix(bed_features_out)
    reports_ch = ann_r1.mix(ann_r3, bed_r1, bed_r2, bed_r3)

    emit:
    features = features_ch
    validation_reports = reports_ch
}
