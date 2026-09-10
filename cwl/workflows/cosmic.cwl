cwlVersion: v1.2
class: Workflow

doc: "COSMIC branch: validate TSV, keep FATHMM-MKL > 0.7, intersect, annotate, summarize, export."

inputs:
  tsv: File
  features: File

outputs:
  validated_tsv:
    type: File
    outputSource: validate_cosmic/validated
  validation_report:
    type: File
    outputSource: validate_cosmic/report
  filtered_bed:
    type: File
    outputSource: validate_filtered_bed/validated_bed
  filtered_tsv:
    type: File
    outputSource: filter_cosmic/tsv
  filter_report:
    type: File
    outputSource: filter_cosmic/report
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
  validate_cosmic:
    run: ../tools/validate_cosmic_tsv.cwl
    in:
      tsv: tsv
    out: [validated, report]

  filter_cosmic:
    run: ../tools/filter_cosmic_fathmm.cwl
    in:
      validated_tsv: validate_cosmic/validated
    out: [bed, tsv, report]

  validate_filtered_bed:
    run: ../tools/validate_bed.cwl
    in:
      bed: filter_cosmic/bed
      min_columns:
        default: 4
      output_name:
        default: cosmic_pathogenic.bed
      report_name:
        default: cosmic_bed_validation_report.tsv
    out: [validated_bed, report]

  intersect:
    run: ../tools/intersect.cwl
    in:
      a: features
      b: validate_filtered_bed/validated_bed
    out: [raw]

  format:
    run: ../tools/format_cosmic_intersections.cwl
    in:
      raw: intersect/raw
    out: [annotations]

  summarize:
    run: ../tools/summarize_annotations.cwl
    in:
      annotations: format/annotations
      summary_name:
        default: gtf_cosmic_summary.tsv
    out: [summary]

  export:
    run: ../tools/export_cosmic_formats.cwl
    in:
      annotations: format/annotations
    out: [sql, ucsc]
