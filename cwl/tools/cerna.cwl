cwlVersion: v1.2
class: CommandLineTool

doc: >
  Run the refactored ceRNA pipeline (parse_tools -> predict_ceRNA ->
  experiment_annotate) from a single JSON manifest. All input paths are
  supplied via the manifest or the --data-root mount.

baseCommand: [python3, /opt/lncbookv3/scripts/cerna_pipeline/pipeline.py]

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
  - cerna_outputs

outputs:
  outputs_dir:
    type: Directory
    outputBinding:
      glob: cerna_outputs
