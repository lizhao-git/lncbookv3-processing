cwlVersion: v1.2
class: CommandLineTool

doc: "Extract gene/transcript/exon/intron features from GTF/GFF3 annotation into a BED-like table."

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/extract_annotation_features.py]

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
  features:
    type:
      type: array
      items: string
    default: [gene, transcript, exon, intron]
    inputBinding:
      prefix: --feature
  output_name:
    type: string
    default: annotation_features.bed

arguments:
  - --output-bed
  - $(inputs.output_name)

outputs:
  features_bed:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
