cwlVersion: v1.2
class: CommandLineTool

doc: >
  Split a VCF into per-contig or per-sample files named
  PREFIX.<key>.vcf[.gz] inside --out-dir. Output compression mirrors the
  input; every output keeps the full meta header.

baseCommand: [python3, /opt/lncbookv3/scripts/format_convert/vcf_split.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest
  InlineJavascriptRequirement: {}

inputs:
  vcf:
    type: File
    inputBinding:
      prefix: --input-vcf
  mode:
    type:
      type: enum
      symbols: [chrom, sample]
    inputBinding:
      prefix: --mode
  prefix:
    type: string
    default: records
    inputBinding:
      prefix: --prefix
  out_dir:
    type: string
    default: split_vcf
  report_name:
    type: string
    default: vcf_split_report.tsv

arguments:
  - --out-dir
  - $(inputs.out_dir)
  - --report
  - $(inputs.report_name)

outputs:
  split_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.out_dir)/*
  report:
    type: File
    outputBinding:
      glob: $(inputs.report_name)
