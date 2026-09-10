cwlVersion: v1.2
class: CommandLineTool

doc: "Summarize feature-level annotations (shared by the ClinVar and COSMIC branches)."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/summarize_annotations.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  annotations:
    type: File
    inputBinding:
      prefix: --input-tsv
  summary_name:
    type: string
    inputBinding:
      prefix: --output-summary

outputs:
  summary:
    type: File
    outputBinding:
      glob: $(inputs.summary_name)
