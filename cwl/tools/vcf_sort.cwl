cwlVersion: v1.2
class: CommandLineTool

doc: >
  Sort VCF records by contig and position. Contig order comes from
  --contig-order (comma-separated), else the .fai index, else the ##contig
  header, else first appearance in the file.

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/vcf_sort.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  vcf:
    type: File
    inputBinding:
      prefix: --input-vcf
  contig_order:
    type: string?
    inputBinding:
      prefix: --contig-order
  fasta_index:
    type: File?
    inputBinding:
      prefix: --fasta-index
  output_name:
    type: string
    default: sorted.vcf
  report_name:
    type: string
    default: vcf_sort_report.tsv

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
