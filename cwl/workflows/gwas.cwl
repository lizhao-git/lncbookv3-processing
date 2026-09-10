cwlVersion: v1.2
class: Workflow

doc: "GWAS Catalog branch: validate TSV, keep p-value < 5e-8, intersect, annotate, summarize, export."

inputs:
  tsv: File
  features: File

outputs:
  validated_tsv:
    type: File
    outputSource: validate_gwas/validated
  validation_report:
    type: File
    outputSource: validate_gwas/report
  filtered_bed:
    type: File
    outputSource: validate_filtered_bed/validated_bed
  filtered_tsv:
    type: File
    outputSource: filter_gwas/tsv
  filter_report:
    type: File
    outputSource: filter_gwas/report
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
  validate_gwas:
    run: ../tools/validate_gwas_tsv.cwl
    in:
      tsv: tsv
    out: [validated, report]

  filter_gwas:
    run: ../tools/filter_gwas_significant.cwl
    in:
      validated_tsv: validate_gwas/validated
    out: [bed, tsv, report]

  validate_filtered_bed:
    run: ../tools/validate_bed.cwl
    in:
      bed: filter_gwas/bed
      min_columns:
        default: 4
      output_name:
        default: gwas_catalog_significant.bed
      report_name:
        default: gwas_bed_validation_report.tsv
    out: [validated_bed, report]

  intersect:
    run: ../tools/intersect.cwl
    in:
      a: features
      b: validate_filtered_bed/validated_bed
    out: [raw]

  format:
    run: ../tools/format_gwas_intersections.cwl
    in:
      raw: intersect/raw
    out: [annotations]

  summarize:
    run: ../tools/summarize_gwas.cwl
    in:
      annotations: format/annotations
    out: [summary]

  export:
    run: ../tools/export_gwas_formats.cwl
    in:
      annotations: format/annotations
    out: [sql, ucsc]
