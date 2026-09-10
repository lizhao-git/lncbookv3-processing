cwlVersion: v1.2
class: CommandLineTool

doc: >
  Summarize VCF content into a metric/value report: record counts,
  multiallelic/no-ALT counts, allele type counts, Ts/Tv ratio (SNVs only),
  and per-contig record counts (contig.<name> rows).

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/vcf_stats.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  vcf:
    type: File
    inputBinding:
      prefix: --input-vcf
  report_name:
    type: string
    default: vcf_stats_report.tsv

arguments:
  - --report
  - $(inputs.report_name)

outputs:
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
