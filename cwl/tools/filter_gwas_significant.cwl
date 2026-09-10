cwlVersion: v1.2
class: CommandLineTool

doc: "Filter GWAS Catalog associations with p-value < 5e-8 into BED/TSV."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/filter_gwas_catalog_significant.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  validated_tsv:
    type: File
    inputBinding:
      prefix: --input-tsv

arguments:
  - --output-bed
  - gwas_catalog_significant.bed
  - --output-tsv
  - gwas_catalog_significant.tsv
  - --report
  - gwas_catalog_filter_report.tsv

outputs:
  bed:
    type: File
    outputBinding:
      glob: gwas_catalog_significant.bed
  tsv:
    type: File
    outputBinding:
      glob: gwas_catalog_significant.tsv
  report:
    type: File
    outputBinding:
      glob: gwas_catalog_filter_report.tsv
