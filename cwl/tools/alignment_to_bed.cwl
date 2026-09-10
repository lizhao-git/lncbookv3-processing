cwlVersion: v1.2
class: CommandLineTool

doc: >
  Convert SAM/BAM alignments into BED6 intervals. BAM input is streamed
  through samtools view; --split-blocks emits one row per aligned block
  instead of spanning introns.

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/alignment_to_bed.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  alignment:
    type: File
    inputBinding:
      prefix: --input-alignment
  min_mapq:
    type: int?
    inputBinding:
      prefix: --min-mapq
  split_blocks:
    type: boolean
    default: false
    inputBinding:
      prefix: --split-blocks
  output_name:
    type: string
    default: alignments.bed
  report_name:
    type: string
    default: alignment_to_bed_report.tsv

arguments:
  - --output-bed
  - $(inputs.output_name)
  - --report
  - $(inputs.report_name)

outputs:
  bed:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
