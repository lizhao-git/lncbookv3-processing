cwlVersion: v1.2
class: CommandLineTool

doc: "Filter COSMIC variants with FATHMM-MKL score > 0.7 into BED/TSV."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/filter_cosmic_fathmm.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  validated_tsv:
    type: File
    inputBinding:
      prefix: --input-tsv

arguments:
  - --output-bed
  - cosmic_pathogenic.bed
  - --output-tsv
  - cosmic_pathogenic.tsv
  - --report
  - cosmic_filter_report.tsv

outputs:
  bed:
    type: File
    outputBinding:
      glob: cosmic_pathogenic.bed
  tsv:
    type: File
    outputBinding:
      glob: cosmic_pathogenic.tsv
  report:
    type: File
    outputBinding:
      glob: cosmic_filter_report.tsv
