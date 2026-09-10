cwlVersion: v1.2
class: CommandLineTool

doc: "Validate the LncBook v3 GTF structure and copy it through as validated.gtf."

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/validate_gtf.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  gtf:
    type: File
    inputBinding:
      prefix: --input-gtf

arguments:
  - --output-gtf
  - validated.gtf
  - --report
  - gtf_validation_report.tsv

outputs:
  validated_gtf:
    type: File
    outputBinding:
      glob: validated.gtf
  report:
    type: File
    outputBinding:
      glob: gtf_validation_report.tsv
