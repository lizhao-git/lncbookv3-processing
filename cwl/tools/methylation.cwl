cwlVersion: v1.2
class: CommandLineTool

doc: >
  Run the refactored methylation pipeline (preprocess -> extract -> difftest ->
  aggregate) from a single JSON manifest. Replaces the legacy
  run_analyze_and_format tool; all paths are supplied via the manifest or the
  --data-root mount (no hard-coded /disk1/... layout).

baseCommand: [python3, /opt/lncbookv3/scripts/methylation/pipeline.py]

requirements:
  DockerRequirement:
    dockerPull: lncbookv3-processing:latest

inputs:
  manifest:
    type: File
    inputBinding:
      prefix: --manifest
  data_root:
    type: Directory?
    inputBinding:
      prefix: --data-root

arguments:
  - --output-dir
  - methylation_outputs
  - --report
  - methylation_report.tsv

outputs:
  outputs_dir:
    type: Directory
    outputBinding:
      glob: methylation_outputs
  report:
    type: File
    outputBinding:
      glob: methylation_report.tsv
