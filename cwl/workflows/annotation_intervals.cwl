cwlVersion: v1.2
class: Workflow

doc: >
  Reusable annotation interval workflow: validate a GTF/GFF3 file, extract
  gene/transcript/exon/intron intervals into a BED-like table, and validate
  the extracted BED intervals.

inputs:
  annotation:
    type: File
  annotation_format:
    type: string
    default: auto
  features:
    type:
      type: array
      items: string
    default: [gene, transcript, exon, intron]

outputs:
  validated_annotation:
    type: File
    outputSource: validate_annotation/validated_annotation
  annotation_validation_report:
    type: File
    outputSource: validate_annotation/report
  features_bed:
    type: File
    outputSource: validate_features/validated_bed
  features_bed_validation_report:
    type: File
    outputSource: validate_features/report

steps:
  validate_annotation:
    run: ../tools/validate_annotation.cwl
    in:
      annotation: annotation
      format: annotation_format
    out: [validated_annotation, report]

  extract_features:
    run: ../tools/extract_annotation_features.cwl
    in:
      annotation: validate_annotation/validated_annotation
      format: annotation_format
      features: features
    out: [features_bed]

  validate_features:
    run: ../tools/validate_bed.cwl
    in:
      bed: extract_features/features_bed
      min_columns:
        default: 10
      output_name:
        default: annotation_features.bed
      report_name:
        default: annotation_features_bed_validation_report.tsv
    out: [validated_bed, report]
