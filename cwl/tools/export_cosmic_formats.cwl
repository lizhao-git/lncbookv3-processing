cwlVersion: v1.2
class: CommandLineTool

doc: "Export COSMIC annotation TSV to SQL and UCSC Genome Browser track BED."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/export_cosmic_annotation_formats.py]

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
  - gtf_cosmic_site_annotations.sql
  - --output-ucsc-bed
  - gtf_cosmic_ucsc_track.bed

outputs:
  sql:
    type: File
    outputBinding:
      glob: gtf_cosmic_site_annotations.sql
  ucsc:
    type: File
    outputBinding:
      glob: gtf_cosmic_ucsc_track.bed
