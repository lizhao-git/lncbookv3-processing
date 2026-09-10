cwlVersion: v1.2
class: CommandLineTool

doc: "Validate and normalize one or more SmProt coordinate files into validated_smprot.tsv."

baseCommand: [python3, /opt/lncbookv3/scripts/smprot_pipeline/smprot_validate.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  smprot_files:
    type: File[]
    inputBinding: {}

arguments:
  - --output-tsv
  - validated_smprot.tsv
  - --report
  - smprot_validation_report.tsv

outputs:
  validated:
    type: File
    outputBinding:
      glob: validated_smprot.tsv
  report:
    type: File
    outputBinding:
      glob: smprot_validation_report.tsv
