cwlVersion: v1.2
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

doc: >
  lncbookv3-processing: map ClinVar / COSMIC / GWAS Catalog / SmProt resources
  onto LncBook v3 lncRNA GTF annotations, with an optional methylation analysis
  branch. Reference branch always runs; optional branches are gated by boolean
  flags using the empty-array scatter idiom.

inputs:
  gtf:
    type: File
  run_clinvar:
    type: boolean
  clinvar_vcf:
    type: File?
  run_cosmic:
    type: boolean
  cosmic_tsv:
    type: File?
  run_gwas_catalog:
    type: boolean
  gwas_tsv:
    type: File?
  run_smprot:
    type: boolean
  smprot_tsv:
    type: File[]
  run_methylation:
    type: boolean
    default: false
  methylation_manifest:
    type: File?
    default: null
  methylation_data_root:
    type: Directory?
    default: null
  run_cerna:
    type: boolean
    default: false
  cerna_manifest:
    type: File?
    default: null
  cerna_data_root:
    type: Directory?
    default: null

outputs:
  # Reference branch (always produced)
  validated_gtf:
    type: File
    outputSource: reference/validated_gtf
  gtf_validation_report:
    type: File
    outputSource: reference/gtf_validation_report
  gtf_features:
    type: File
    outputSource: reference/features
  gtf_features_validation_report:
    type: File
    outputSource: reference/features_validation_report

  # ClinVar branch
  clinvar_validated_vcf:
    type: File[]
    outputSource: clinvar/validated_vcf
  clinvar_validation_report:
    type: File[]
    outputSource: clinvar/validation_report
  clinvar_filtered_bed:
    type: File[]
    outputSource: clinvar/filtered_bed
  clinvar_filtered_vcf:
    type: File[]
    outputSource: clinvar/filtered_vcf
  clinvar_filter_report:
    type: File[]
    outputSource: clinvar/filter_report
  clinvar_raw_intersections:
    type: File[]
    outputSource: clinvar/raw_intersections
  clinvar_annotations:
    type: File[]
    outputSource: clinvar/annotations
  clinvar_summary:
    type: File[]
    outputSource: clinvar/summary
  clinvar_sql:
    type: File[]
    outputSource: clinvar/sql
  clinvar_ucsc:
    type: File[]
    outputSource: clinvar/ucsc

  # COSMIC branch
  cosmic_validated_tsv:
    type: File[]
    outputSource: cosmic/validated_tsv
  cosmic_validation_report:
    type: File[]
    outputSource: cosmic/validation_report
  cosmic_filtered_bed:
    type: File[]
    outputSource: cosmic/filtered_bed
  cosmic_filtered_tsv:
    type: File[]
    outputSource: cosmic/filtered_tsv
  cosmic_filter_report:
    type: File[]
    outputSource: cosmic/filter_report
  cosmic_raw_intersections:
    type: File[]
    outputSource: cosmic/raw_intersections
  cosmic_annotations:
    type: File[]
    outputSource: cosmic/annotations
  cosmic_summary:
    type: File[]
    outputSource: cosmic/summary
  cosmic_sql:
    type: File[]
    outputSource: cosmic/sql
  cosmic_ucsc:
    type: File[]
    outputSource: cosmic/ucsc

  # GWAS Catalog branch
  gwas_validated_tsv:
    type: File[]
    outputSource: gwas/validated_tsv
  gwas_validation_report:
    type: File[]
    outputSource: gwas/validation_report
  gwas_filtered_bed:
    type: File[]
    outputSource: gwas/filtered_bed
  gwas_filtered_tsv:
    type: File[]
    outputSource: gwas/filtered_tsv
  gwas_filter_report:
    type: File[]
    outputSource: gwas/filter_report
  gwas_raw_intersections:
    type: File[]
    outputSource: gwas/raw_intersections
  gwas_annotations:
    type: File[]
    outputSource: gwas/annotations
  gwas_summary:
    type: File[]
    outputSource: gwas/summary
  gwas_sql:
    type: File[]
    outputSource: gwas/sql
  gwas_ucsc:
    type: File[]
    outputSource: gwas/ucsc

  # SmProt branch
  smprot_validated_tsv:
    type: File[]
    outputSource: smprot/validated_tsv
  smprot_validation_report:
    type: File[]
    outputSource: smprot/validation_report
  smprot_filtered_bed:
    type: File[]
    outputSource: smprot/filtered_bed
  smprot_filtered_tsv:
    type: File[]
    outputSource: smprot/filtered_tsv
  smprot_filter_report:
    type: File[]
    outputSource: smprot/filter_report
  smprot_raw_intersections:
    type: File[]
    outputSource: smprot/raw_intersections
  smprot_annotations:
    type: File[]
    outputSource: smprot/annotations
  smprot_mapping_report:
    type: File[]
    outputSource: smprot/mapping_report
  smprot_summary:
    type: File[]
    outputSource: smprot/summary
  smprot_sql:
    type: File[]
    outputSource: smprot/sql
  smprot_ucsc:
    type: File[]
    outputSource: smprot/ucsc

  # Methylation analysis branch
  methylation_outputs:
    type: Directory[]
    outputSource: methylation/outputs_dir
  methylation_report:
    type: File[]
    outputSource: methylation/report

  # ceRNA analysis branch
  cerna_outputs:
    type: Directory[]
    outputSource: cerna/outputs_dir

steps:
  reference:
    run: workflows/reference.cwl
    in:
      gtf: gtf
    out: [validated_gtf, gtf_validation_report, features, features_validation_report]

  clinvar_gate:
    run: tools/when_file.cwl
    in:
      flag: run_clinvar
      value: clinvar_vcf
    out: [result]

  clinvar:
    run: workflows/clinvar.cwl
    in:
      vcf: clinvar_gate/result
      features: reference/features
    scatter: vcf
    out: [validated_vcf, validation_report, filtered_bed, filtered_vcf,
          filter_report, raw_intersections, annotations, summary, sql, ucsc]

  cosmic_gate:
    run: tools/when_file.cwl
    in:
      flag: run_cosmic
      value: cosmic_tsv
    out: [result]

  cosmic:
    run: workflows/cosmic.cwl
    in:
      tsv: cosmic_gate/result
      features: reference/features
    scatter: tsv
    out: [validated_tsv, validation_report, filtered_bed, filtered_tsv,
          filter_report, raw_intersections, annotations, summary, sql, ucsc]

  gwas_gate:
    run: tools/when_file.cwl
    in:
      flag: run_gwas_catalog
      value: gwas_tsv
    out: [result]

  gwas:
    run: workflows/gwas.cwl
    in:
      tsv: gwas_gate/result
      features: reference/features
    scatter: tsv
    out: [validated_tsv, validation_report, filtered_bed, filtered_tsv,
          filter_report, raw_intersections, annotations, summary, sql, ucsc]

  smprot_gate:
    run: tools/when_files.cwl
    in:
      flag: run_smprot
      files: smprot_tsv
    out: [result]

  smprot:
    run: workflows/smprot.cwl
    in:
      files: smprot_tsv
      features: reference/features
      marker: smprot_gate/result
    scatter: marker
    out: [validated_tsv, validation_report, filtered_bed, filtered_tsv,
          filter_report, raw_intersections, annotations, mapping_report,
          summary, sql, ucsc]

  methylation_gate:
    run: tools/when_file.cwl
    in:
      flag: run_methylation
      value: methylation_manifest
    out: [result]

  methylation:
    run: workflows/methylation.cwl
    in:
      manifest: methylation_gate/result
      data_root: methylation_data_root
    scatter: manifest
    out: [outputs_dir, report]

  cerna_gate:
    run: tools/when_file.cwl
    in:
      flag: run_cerna
      value: cerna_manifest
    out: [result]

  cerna:
    run: workflows/cerna.cwl
    in:
      manifest: cerna_gate/result
      data_root: cerna_data_root
    scatter: manifest
    out: [outputs_dir]
