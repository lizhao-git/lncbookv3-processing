cwlVersion: v1.2
class: CommandLineTool

doc: "Filter ClinVar variants with definite clinical labels into BED/VCF."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/filter_clinvar_labels.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  validated_vcf:
    type: File
    inputBinding:
      prefix: --input-vcf

arguments:
  - --output-bed
  - clinvar_definite_labels.bed
  - --output-vcf
  - clinvar_definite_labels.vcf
  - --report
  - clinvar_filter_report.tsv

outputs:
  bed:
    type: File
    outputBinding:
      glob: clinvar_definite_labels.bed
  vcf:
    type: File
    outputBinding:
      glob: clinvar_definite_labels.vcf
  report:
    type: File
    outputBinding:
      glob: clinvar_filter_report.tsv
