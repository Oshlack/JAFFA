#!/usr/bin/env cwl-runner
#
# Author: Rebecca Evans (rebecca.louise.evans@gmail.com)

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "bpipe", "run" ]
label: JAFFA
doc: "Jaffa: Just Another Fusion Finding Algorithm"

stdout: log.txt

# It is preferred to have the docker requirement as a hint
hints:
#  InitialWorkDirRequirement:
#    listing:
  DockerRequirement:
    dockerPull: beccyl/jaffa:1.09_dev
    dockerOutputDirectory: /usr/local/jaffa/batch

requirements:
  - class: InlineJavascriptRequirement
#  - class: ResourceRequirement
#    coresMin: 8
#    ramMin: 60000

inputs:

  threads:
    type: int?
    doc: | 
      Change number of threads used
    inputBinding:
      position: 2
      prefix: -n

  refBase:
    type: Directory
    inputBinding:
      position: 4
      prefix: refBase=
      separate: false

  genomeFasta:
    type: File
    inputBinding:
      position: 6
      prefix: genomeFasta=
      separate: false

  fastqInputFormat:
    type: string
    inputBinding:
      position: 8
      prefix: fastqInputFormat=
      separate: false

  genome:
    type: string
    inputBinding:
      position: 10
      prefix: genome=
      separate: false

  annotation:
    type: string
    inputBinding:
      position: 12
      prefix: annotation=
      separate: false

## path to groovy script within docker image ie. /opt/JAFFA/JAFFA_direct.groovy
  pipeline:
    type: string
    inputBinding:
      position: 13

  fastq1:
    type: File
    inputBinding:
      position: 14

  fastq2:
    type: File
    inputBinding:
      position: 15

  
outputs:
  output:
    type: stdout
  
  jaffa_fusions_fasta:
    type: File
    outputBinding:
      glob: "jaffa_results.fasta"

  jaffa_fusions_list:
    type: File
    outputBinding:
      glob: "jaffa_results.csv"

arguments:
  - valueFrom: -p
    position: 3
  - valueFrom: -p
    position: 5
  - valueFrom: -p
    position: 7
  - valueFrom: -p
    position: 9
  - valueFrom: -p
    position: 11
