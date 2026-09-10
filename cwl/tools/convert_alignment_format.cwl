cwlVersion: v1.2
class: CommandLineTool

doc: >
  Convert between SAM and BAM with samtools. The target format is inferred
  from the output file extension unless --to is set; sorting and BAM
  indexing are optional. Requires samtools (shipped in the image).

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/convert_alignment_format.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  alignment:
    type: File
    inputBinding:
      prefix: --input-alignment
  to_format:
    type:
      type: enum
      symbols: [auto, sam, bam]
    default: auto
    inputBinding:
      prefix: --to
  sort:
    type: boolean
    default: false
    inputBinding:
      prefix: --sort
  index:
    type: boolean
    default: false
    inputBinding:
      prefix: --index
  threads:
    type: int
    default: 0
    inputBinding:
      prefix: --threads
  output_name:
    type: string
    default: converted_alignment.bam

arguments:
  - --output-alignment
  - $(inputs.output_name)

outputs:
  alignment:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
