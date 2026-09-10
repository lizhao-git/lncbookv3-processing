cwlVersion: v1.2
class: CommandLineTool

doc: "Convert a BED interval file into GFF3 (0-based BED to 1-based GFF3)."

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/bed_to_gff3.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  bed:
    type: File
    inputBinding:
      prefix: --input-bed
  feature:
    type: string
    default: region
    inputBinding:
      prefix: --feature
  source:
    type: string
    default: "."
    inputBinding:
      prefix: --source
  output_name:
    type: string
    default: bed_annotation.gff3

arguments:
  - --output-gff3
  - $(inputs.output_name)

outputs:
  gff3:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
