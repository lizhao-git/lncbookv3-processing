cwlVersion: v1.2
class: CommandLineTool

doc: "Validate a BED coordinate file and copy it through."

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/validate_bed.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  bed:
    type: File
    inputBinding:
      prefix: --input-bed
  min_columns:
    type: int
    default: 3
    inputBinding:
      prefix: --min-columns
  output_name:
    type: string
    default: validated.bed
  report_name:
    type: string
    default: bed_validation_report.tsv

arguments:
  - --output-bed
  - $(inputs.output_name)
  - --report
  - $(inputs.report_name)

outputs:
  validated_bed:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
