cwlVersion: v1.2
class: CommandLineTool

doc: >
  Convert a BED file into a bigBed track. The bedN type is inferred from
  the uniform column count unless --bed-type is set; rows are sorted into
  the chrom.sizes order and converted with bedToBigBed. Requires the UCSC
  kent tools (shipped in the image).

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/bed_to_bigbed.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  bed:
    type: File
    inputBinding:
      prefix: --input-bed
  chrom_sizes:
    type: File
    inputBinding:
      prefix: --chrom-sizes
  bed_type:
    type: string?
    inputBinding:
      prefix: --bed-type
  output_name:
    type: string
    default: track.bigBed
  report_name:
    type: string
    default: bed_to_bigbed_report.tsv

arguments:
  - --output-bigbed
  - $(inputs.output_name)
  - --report
  - $(inputs.report_name)

outputs:
  bigbed:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
