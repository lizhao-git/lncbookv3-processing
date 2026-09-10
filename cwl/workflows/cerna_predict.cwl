cwlVersion: v1.2
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}

doc: >
  ceRNA target-prediction workflow: prepare lncRNA/mature-miRNA sequences from
  the GTF + reference genome + miRBase mature.fa, then run miRanda / RNAhybrid /
  TargetScan and emit the normalized prediction files for the downstream ceRNA
  pipeline.

inputs:
  manifest:
    type: File
  data_root:
    type: Directory?

outputs:
  outputs_dir:
    type: Directory
    outputSource: predict/outputs_dir

steps:
  predict:
    run: ../tools/cerna_predict.cwl
    in:
      manifest: manifest
      data_root: data_root
    out: [outputs_dir]
