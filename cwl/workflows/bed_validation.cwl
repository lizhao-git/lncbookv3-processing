cwlVersion: v1.2
class: Workflow

doc: "Reusable BED validation workflow for 0-based half-open interval files."

inputs:
  bed:
    type: File
  min_columns:
    type: int
    default: 3

outputs:
  validated_bed:
    type: File
    outputSource: validate_bed/validated_bed
  bed_validation_report:
    type: File
    outputSource: validate_bed/report

steps:
  validate_bed:
    run: ../tools/validate_bed.cwl
    in:
      bed: bed
      min_columns: min_columns
    out: [validated_bed, report]
