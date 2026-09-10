cwlVersion: v1.2
class: CommandLineTool

doc: "Summarize retained SmProt-lncRNA mappings."

baseCommand: [python3, /opt/lncbookv3/scripts/smprot_pipeline/summarize_smprot_annotations.py]

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
  - gtf_smprot_summary.tsv

outputs:
  summary:
    type: File
    outputBinding:
      glob: gtf_smprot_summary.tsv
