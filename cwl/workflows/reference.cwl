cwlVersion: v1.2
class: Workflow

doc: "Reference branch: validate the LncBook annotation and extract reusable BED features."

inputs:
  gtf: File
  annotation_format:
    type: string
    default: gtf

outputs:
  validated_gtf:
    type: File
    outputSource: validate_annotation/validated_annotation
  gtf_validation_report:
    type: File
    outputSource: validate_annotation/report
  features:
    type: File
    outputSource: validate_features/validated_bed
  features_validation_report:
    type: File
    outputSource: validate_features/report

steps:
  validate_annotation:
    run: ../tools/validate_annotation.cwl
    in:
      annotation: gtf
      format: annotation_format
      output_name:
        default: validated.gtf
      report_name:
        default: gtf_validation_report.tsv
    out: [validated_annotation, report]

  extract_features:
    run: ../tools/extract_annotation_features.cwl
    in:
      annotation: validate_annotation/validated_annotation
      format: annotation_format
      output_name:
        default: gtf_features.bed
    out: [features_bed]

  validate_features:
    run: ../tools/validate_bed.cwl
    in:
      bed: extract_features/features_bed
      min_columns:
        default: 10
      output_name:
        default: gtf_features.bed
      report_name:
        default: gtf_features_validation_report.tsv
    out: [validated_bed, report]
