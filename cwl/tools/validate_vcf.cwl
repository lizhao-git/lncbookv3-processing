cwlVersion: v1.2
class: CommandLineTool

doc: "Validate a VCF file (headers, record structure, REF/ALT/QUAL/FILTER/INFO, FORMAT/samples) and copy it through as validated_vcf.txt."

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/validate_vcf.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  vcf:
    type: File
    inputBinding:
      prefix: --input-vcf
  output_name:
    type: string
    default: validated_vcf.txt
  report_name:
    type: string
    default: vcf_validation_report.tsv

arguments:
  - --output-vcf
  - $(inputs.output_name)
  - --report
  - $(inputs.report_name)

outputs:
  validated_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
