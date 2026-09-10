cwlVersion: v1.2
class: CommandLineTool

doc: >
  Convert a bedGraph file into a bigWig track. Rows are sorted into the
  chrom.sizes order and converted with bedGraphToBigWig; --clip clamps
  regions that overrun chromosome bounds. Requires the UCSC kent tools
  (shipped in the image).

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/bedgraph_to_bigwig.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  bedgraph:
    type: File
    inputBinding:
      prefix: --input-bedgraph
  chrom_sizes:
    type: File
    inputBinding:
      prefix: --chrom-sizes
  clip:
    type: boolean
    default: false
    inputBinding:
      prefix: --clip
  output_name:
    type: string
    default: track.bigWig
  report_name:
    type: string
    default: bedgraph_to_bigwig_report.tsv

arguments:
  - --output-bigwig
  - $(inputs.output_name)
  - --report
  - $(inputs.report_name)

outputs:
  bigwig:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
