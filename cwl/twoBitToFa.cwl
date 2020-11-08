#!/usr/bin/env cwl-runner
#
# Author: Rebecca Evans (rebecca.louise.evans@gmail.com)

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [twoBitToFa]
#stdout: $(inputs.output)

doc: "Convert jaffa fusion output to bedpe format"

hints:
  DockerRequirement:
    dockerPull: beccyl/jaffa:1.09_dev

requirements:
  - class: InlineJavascriptRequirement

inputs:

  input:
    type: File
    inputBinding:
      position: 1

  output:
    type: string
    inputBinding:
      position: 2

outputs:

  output:
    type: File
    outputBinding:
      glob: $(inputs.output)

