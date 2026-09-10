cwlVersion: v1.2
class: CommandLineTool

doc: "Validate GTF/GFF3 annotation and copy it through as validated_annotation.txt."

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/validate_annotation.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  annotation:
    type: File
    inputBinding:
      prefix: --input-annotation
  format:
    type: string
    default: auto
    inputBinding:
      prefix: --format
  output_name:
    type: string
    default: validated_annotation.txt
  report_name:
    type: string
    default: annotation_validation_report.tsv

arguments:
  - --output-annotation
  - $(inputs.output_name)
  - --report
  - $(inputs.report_name)

outputs:
  validated_annotation:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
