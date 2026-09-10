cwlVersion: v1.2
class: Workflow

doc: >
  SmProt branch: validate/normalize coordinate files, emit BED, intersect,
  keep unique in-transcript mappings, summarize, export. The `marker` input is
  only used as a scatter trigger by the top-level workflow (zero or one runs)
  and is otherwise ignored.

inputs:
  files: File[]
  features: File
  marker: File

outputs:
  validated_tsv:
    type: File
    outputSource: validate_smprot/validated
  validation_report:
    type: File
    outputSource: validate_smprot/report
  filtered_bed:
    type: File
    outputSource: validate_filtered_bed/validated_bed
  filtered_tsv:
    type: File
    outputSource: filter_smprot/tsv
  filter_report:
    type: File
    outputSource: filter_smprot/report
  raw_intersections:
    type: File
    outputSource: intersect/raw
  annotations:
    type: File
    outputSource: format/annotations
  mapping_report:
    type: File
    outputSource: format/report
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
  validate_smprot:
    run: ../tools/validate_smprot_tsv.cwl
    in:
      smprot_files: files
    out: [validated, report]

  filter_smprot:
    run: ../tools/filter_smprot_records.cwl
    in:
      validated_tsv: validate_smprot/validated
    out: [bed, tsv, report]

  validate_filtered_bed:
    run: ../tools/validate_bed.cwl
    in:
      bed: filter_smprot/bed
      min_columns:
        default: 4
      output_name:
        default: smprot_filtered.bed
      report_name:
        default: smprot_bed_validation_report.tsv
    out: [validated_bed, report]

  intersect:
    run: ../tools/intersect.cwl
    in:
      a: features
      b: validate_filtered_bed/validated_bed
    out: [raw]

  format:
    run: ../tools/format_smprot_intersections.cwl
    in:
      raw: intersect/raw
    out: [annotations, report]

  summarize:
    run: ../tools/summarize_smprot.cwl
    in:
      annotations: format/annotations
    out: [summary]

  export:
    run: ../tools/export_smprot_formats.cwl
    in:
      annotations: format/annotations
    out: [sql, ucsc]
