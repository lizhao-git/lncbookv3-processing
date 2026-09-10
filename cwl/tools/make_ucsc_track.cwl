cwlVersion: v1.2
class: CommandLineTool

doc: >
  Add or replace a UCSC track line on a BED/bedGraph file. Existing track
  lines are replaced and browser lines dropped, so the output carries
  exactly one track definition.

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/make_ucsc_track.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  bed:
    type: File
    inputBinding:
      prefix: --input-bed
  track_name:
    type: string
    inputBinding:
      prefix: --track-name
  description:
    type: string?
    inputBinding:
      prefix: --description
  track_type:
    type: string?
    inputBinding:
      prefix: --track-type
  color:
    type: string?
    inputBinding:
      prefix: --color
  visibility:
    type: string?
    inputBinding:
      prefix: --visibility
  priority:
    type: int?
    inputBinding:
      prefix: --priority
  output_name:
    type: string
    default: track.bed
  report_name:
    type: string
    default: ucsc_track_report.tsv

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
