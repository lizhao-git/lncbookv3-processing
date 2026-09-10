cwlVersion: v1.2
class: Workflow

doc: >
  ceRNA analysis branch (refactored Python 3 pipeline). Takes a JSON manifest
  plus an optional data-root Directory that anchors the manifest's relative
  paths, and produces predicted ceRNA interactions, experimental mappings and
  disease annotations.

inputs:
  manifest: File
  data_root: Directory?

outputs:
  outputs_dir:
    type: Directory
    outputSource: run_cerna/outputs_dir

steps:
  run_cerna:
    run: ../tools/cerna.cwl
    in:
      manifest: manifest
      data_root: data_root
    out: [outputs_dir]
