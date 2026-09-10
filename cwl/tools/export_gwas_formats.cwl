cwlVersion: v1.2
class: CommandLineTool

doc: "Export GWAS Catalog annotation TSV to SQL and UCSC Genome Browser track BED."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/export_gwas_catalog_annotation_formats.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  annotations:
    type: File
    inputBinding:
      prefix: --input-tsv

arguments:
  - --output-sql
  - gtf_gwas_catalog_site_annotations.sql
  - --output-ucsc-bed
  - gtf_gwas_catalog_ucsc_track.bed

outputs:
  sql:
    type: File
    outputBinding:
      glob: gtf_gwas_catalog_site_annotations.sql
  ucsc:
    type: File
    outputBinding:
      glob: gtf_gwas_catalog_ucsc_track.bed
