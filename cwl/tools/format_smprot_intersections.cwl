cwlVersion: v1.2
class: CommandLineTool

doc: "Retain SmProt mappings entirely and uniquely within lncRNA transcripts."

baseCommand: [python3, /opt/lncbookv3/scripts/smprot_pipeline/format_smprot_intersections.py]

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
  - gtf_smprot_site_annotations.tsv
  - --report
  - smprot_mapping_report.tsv

outputs:
  annotations:
    type: File
    outputBinding:
      glob: gtf_smprot_site_annotations.tsv
  report:
    type: File
    outputBinding:
      glob: smprot_mapping_report.tsv
