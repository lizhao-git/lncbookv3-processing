cwlVersion: v1.2
class: Workflow

doc: >
  Methylation analysis branch (refactored Python 3 pipeline). Takes a JSON
  manifest plus an optional data-root Directory that anchors the manifest's
  relative paths, and produces the aggregate significance matrix, labels and
  the database-import table.

inputs:
  manifest: File
  data_root: Directory?

outputs:
  outputs_dir:
    type: Directory
    outputSource: run_methylation/outputs_dir
  report:
    type: File
    outputSource: run_methylation/report

steps:
  run_methylation:
    run: ../tools/methylation.cwl
    in:
      manifest: manifest
      data_root: data_root
    out: [outputs_dir, report]
