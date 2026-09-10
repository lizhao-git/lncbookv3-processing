cwlVersion: v1.2
class: CommandLineTool

doc: "Summarize feature-level GWAS Catalog annotations."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/summarize_gwas_catalog_annotations.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  annotations:
    type: File
    inputBinding:
      prefix: --input-tsv

arguments:
  - --output-summary
  - gtf_gwas_catalog_summary.tsv

outputs:
  summary:
    type: File
    outputBinding:
      glob: gtf_gwas_catalog_summary.tsv
