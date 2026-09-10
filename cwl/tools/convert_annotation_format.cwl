cwlVersion: v1.2
class: CommandLineTool

doc: >
  Convert between GTF and GFF3 annotation formats. GTF -> GFF3 adds
  ID/Parent attributes for the gene/transcript/child hierarchy; GFF3 -> GTF
  restores gene_id/transcript_id from the Parent chain.

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/convert_annotation_format.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  annotation:
    type: File
    inputBinding:
      prefix: --input-annotation
  from_format:
    type: string
    default: auto
    inputBinding:
      prefix: --from
  to_format:
    type:
      type: enum
      symbols: [gtf, gff3]
    inputBinding:
      prefix: --to
  output_name:
    type: string
    default: '$(inputs.to_format == "gtf" ? "converted_annotation.gtf" : "converted_annotation.gff3")'
  report_name:
    type: string
    default: annotation_conversion_report.tsv

arguments:
  - --output-annotation
  - $(inputs.output_name)
  - --report
  - $(inputs.report_name)

outputs:
  converted_annotation:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
