cwlVersion: v1.2
class: CommandLineTool

doc: >
  Filter VCF records by region, QUAL, INFO/DP, allele type, or FILTER
  status. Checks run in a fixed order (region, QUAL, DP, type, FILTER) and
  drop a record at the first failing check; the report counts removals per
  reason.

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/vcf_filter.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  vcf:
    type: File
    inputBinding:
      prefix: --input-vcf
  regions:
    type: string[]?
    inputBinding:
      prefix: --region
  min_qual:
    type: float?
    inputBinding:
      prefix: --min-qual
  min_dp:
    type: int?
    inputBinding:
      prefix: --min-dp
  keep_types:
    type:
      type: array
      items:
        type: enum
        symbols: [snv, mnv, indel, other]
    inputBinding:
      prefix: --keep-type
  require_pass:
    type: boolean
    default: false
    inputBinding:
      prefix: --require-pass
  output_name:
    type: string
    default: filtered.vcf
  report_name:
    type: string
    default: vcf_filter_report.tsv

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
