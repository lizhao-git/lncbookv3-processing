cwlVersion: v1.2
class: CommandLineTool

doc: >
  Run the ceRNA target-prediction stage (prepare_sequences -> predict_targets)
  from a single JSON manifest. Produces the three normalized prediction files
  (miranda / rnahybrid / targetscan) consumed by the downstream ceRNA pipeline.

baseCommand: [python3, /opt/lncbookv3/scripts/cerna_pipeline/pipeline_predict.py]

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
  - cerna_predict

outputs:
  outputs_dir:
    type: Directory
    outputBinding:
      glob: cerna_predict
