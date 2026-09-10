cwlVersion: v1.2
class: Workflow

doc: "ClinVar branch: validate VCF, keep definite labels, intersect, annotate, summarize, export."

inputs:
  vcf: File
  features: File

outputs:
  validated_vcf:
    type: File
    outputSource: validate_clinvar/validated
  validation_report:
    type: File
    outputSource: validate_clinvar/report
  filtered_bed:
    type: File
    outputSource: validate_filtered_bed/validated_bed
  filtered_vcf:
    type: File
    outputSource: filter_clinvar/vcf
  filter_report:
    type: File
    outputSource: filter_clinvar/report
  raw_intersections:
    type: File
    outputSource: intersect/raw
  annotations:
    type: File
    outputSource: format/annotations
  summary:
    type: File
    outputSource: summarize/summary
  sql:
    type: File
    outputSource: export/sql
  ucsc:
    type: File
    outputSource: export/ucsc

steps:
  validate_clinvar:
    run: ../tools/validate_clinvar_vcf.cwl
    in:
      vcf: vcf
    out: [validated, report]

  filter_clinvar:
    run: ../tools/filter_clinvar_labels.cwl
    in:
      validated_vcf: validate_clinvar/validated
    out: [bed, vcf, report]

  validate_filtered_bed:
    run: ../tools/validate_bed.cwl
    in:
      bed: filter_clinvar/bed
      min_columns:
        default: 4
      output_name:
        default: clinvar_definite_labels.bed
      report_name:
        default: clinvar_bed_validation_report.tsv
    out: [validated_bed, report]

  intersect:
    run: ../tools/intersect.cwl
    in:
      a: features
      b: validate_filtered_bed/validated_bed
    out: [raw]

  format:
    run: ../tools/format_clinvar_intersections.cwl
    in:
      raw: intersect/raw
    out: [annotations]

  summarize:
    run: ../tools/summarize_annotations.cwl
    in:
      annotations: format/annotations
      summary_name:
        default: gtf_clinvar_summary.tsv
    out: [summary]

  export:
    run: ../tools/export_clinvar_formats.cwl
    in:
      annotations: format/annotations
    out: [sql, ucsc]
