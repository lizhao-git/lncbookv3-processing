cwlVersion: v1.2
class: CommandLineTool

doc: "Filter/normalize SmProt coordinates and emit a BED file."

baseCommand: [python3, /opt/lncbookv3/scripts/smprot_pipeline/filter_smprot_records.py]

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
  - smprot_filtered.bed
  - --output-tsv
  - smprot_filtered.tsv
  - --report
  - smprot_filter_report.tsv

outputs:
  bed:
    type: File
    outputBinding:
      glob: smprot_filtered.bed
  tsv:
    type: File
    outputBinding:
      glob: smprot_filtered.tsv
  report:
    type: File
    outputBinding:
      glob: smprot_filter_report.tsv
