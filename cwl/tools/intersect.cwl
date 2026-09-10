cwlVersion: v1.2
class: CommandLineTool

doc: "Intersect two BED files with `bedtools intersect -wa -wb` and capture stdout."

baseCommand: bedtools

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  a:
    type: File
    inputBinding:
      prefix: -a
  b:
    type: File
    inputBinding:
      prefix: -b

arguments:
  - intersect
  - -wa
  - -wb

stdout: raw_intersections.tsv

outputs:
  raw:
    type: File
    outputBinding:
      glob: raw_intersections.tsv
