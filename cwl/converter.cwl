#!/usr/bin/env cwl-runner
# Author: Rebecca Evans (rebecca.louise.evans@gmail.com)

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [convert_jaffa_to_bedpe.py]
stdout: $(inputs.output)

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

  output: string

outputs:

  fusionout:
    type: File
    outputBinding:
      glob: $(inputs.output)

