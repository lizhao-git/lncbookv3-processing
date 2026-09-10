cwlVersion: v1.2
class: ExpressionTool

doc: >
  Gate a single-file branch on a boolean flag. Returns a one-element array when
  the flag is true, otherwise an empty array. Scattering a sub-workflow over the
  result runs it zero or one times (the standard CWL idiom for optional steps).

requirements:
  InlineJavascriptRequirement: {}

inputs:
  flag:
    type: boolean
  value:
    type: File?

outputs:
  result:
    type: File[]

expression: |
  ${
    if (inputs.flag && inputs.value) {
      return {result: [inputs.value]};
    }
    return {result: []};
  }
