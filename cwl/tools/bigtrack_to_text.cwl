cwlVersion: v1.2
class: CommandLineTool

doc: >
  Convert a bigBed track back to BED or a bigWig track back to bedGraph.
  The tool is chosen from the input extension (.bigBed/.bb vs
  .bigWig/.bw). Requires the UCSC kent tools (shipped in the image).

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/bigtrack_to_text.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  track:
    type: File
    inputBinding:
      prefix: --input-track
  output_name:
    type: string
    default: track.bed
  report_name:
    type: string
    default: bigtrack_to_text_report.tsv

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
