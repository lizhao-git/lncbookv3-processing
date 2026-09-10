cwlVersion: v1.2
class: ExpressionTool

doc: >
  Gate a multi-file branch on a boolean flag. Returns a one-element marker array
  when the flag is true (the first file is used only as a scatter marker), or an
  empty array otherwise. Scattering the branch workflow over the marker runs it
  zero or one times while still passing the full file list as a constant input.

requirements:
  InlineJavascriptRequirement: {}

inputs:
  flag:
    type: boolean
  files:
    type: File[]?

outputs:
  result:
    type: File[]

expression: |
  ${
    if (inputs.flag && inputs.files && inputs.files.length > 0) {
      return {result: [inputs.files[0]]};
    }
    return {result: []};
  }
