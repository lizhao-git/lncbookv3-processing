cwlVersion: v1.2
class: CommandLineTool

doc: "Validate the GWAS Catalog TSV structure and copy it through as validated_gwas_catalog.tsv."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/validate_gwas_catalog_tsv.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  tsv:
    type: File
    inputBinding:
      prefix: --input-tsv

arguments:
  - --output-tsv
  - validated_gwas_catalog.tsv
  - --report
  - gwas_catalog_validation_report.tsv

outputs:
  validated:
    type: File
    outputBinding:
      glob: validated_gwas_catalog.tsv
  report:
    type: File
    outputBinding:
      glob: gwas_catalog_validation_report.tsv
