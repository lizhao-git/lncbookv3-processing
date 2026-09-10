cwlVersion: v1.2
class: CommandLineTool

doc: "Format ClinVar bedtools intersections into an annotation TSV."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/format_intersections.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  raw:
    type: File
    inputBinding:
      prefix: --input-intersections

arguments:
  - --output-tsv
  - gtf_clinvar_site_annotations.tsv

outputs:
  annotations:
    type: File
    outputBinding:
      glob: gtf_clinvar_site_annotations.tsv
