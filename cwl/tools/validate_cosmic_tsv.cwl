cwlVersion: v1.2
class: CommandLineTool

doc: "Validate the COSMIC TSV structure and copy it through as validated_cosmic.tsv."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/validate_cosmic_tsv.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  tsv:
    type: File
    inputBinding:
      prefix: --input-tsv

arguments:
  - --output-tsv
  - validated_cosmic.tsv
  - --report
  - cosmic_validation_report.tsv

outputs:
  validated:
    type: File
    outputBinding:
      glob: validated_cosmic.tsv
  report:
    type: File
    outputBinding:
      glob: cosmic_validation_report.tsv
