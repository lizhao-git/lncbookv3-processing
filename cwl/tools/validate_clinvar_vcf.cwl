cwlVersion: v1.2
class: CommandLineTool

doc: "Validate the ClinVar VCF structure and copy it through as validated_clinvar.vcf."

baseCommand: [python3, /opt/lncbookv3/scripts/variant_pipeline/validate_clinvar_vcf.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  vcf:
    type: File
    inputBinding:
      prefix: --input-vcf

arguments:
  - --output-vcf
  - validated_clinvar.vcf
  - --report
  - clinvar_vcf_validation_report.tsv

outputs:
  validated:
    type: File
    outputBinding:
      glob: validated_clinvar.vcf
  report:
    type: File
    outputBinding:
      glob: clinvar_vcf_validation_report.tsv
