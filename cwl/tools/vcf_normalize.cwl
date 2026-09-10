cwlVersion: v1.2
class: CommandLineTool

doc: >
  Normalize a VCF by splitting multiallelic records into biallelic ones.
  GT is recoded per ALT, INFO/FORMAT values are subset by their declared
  Number (A/R/G) with count-based fallback, and the report counts
  multiallelic splits plus unresolved genotype fields.

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/vcf_normalize.py]

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
    default: normalized.vcf
  report_name:
    type: string
    default: vcf_normalize_report.tsv

arguments:
  - --output-vcf
  - $(inputs.output_name)
  - --report
  - $(inputs.report_name)

outputs:
  vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
