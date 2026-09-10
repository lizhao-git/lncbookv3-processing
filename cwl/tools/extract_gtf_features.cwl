cwlVersion: v1.2
class: CommandLineTool

doc: "Extract gene/transcript/exon/intron features from the GTF into a BED file."

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/extract_gtf_features.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  validated_gtf:
    type: File
    inputBinding:
      prefix: --input-gtf

arguments:
  - --output-bed
  - gtf_features.bed

outputs:
  features:
    type: File
    outputBinding:
      glob: gtf_features.bed
