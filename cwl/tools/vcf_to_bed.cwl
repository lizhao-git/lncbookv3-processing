cwlVersion: v1.2
class: CommandLineTool

doc: "Convert a VCF into BED intervals (one row per record, ALT kept as last column)."

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/vcf_to_bed.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  vcf:
    type: File
    inputBinding:
      prefix: --input-vcf
  name_field:
    type:
      type: enum
      symbols: [id, chrom_pos]
    default: id
    inputBinding:
      prefix: --name-field
  output_name:
    type: string
    default: variants.bed

arguments:
  - --output-bed
  - $(inputs.output_name)

outputs:
  bed:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
