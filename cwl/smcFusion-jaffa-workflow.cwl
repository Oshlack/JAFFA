#!/usr/bin/env cwl-runner
#
# Author: Rebecca Evans (rebecca.louise.evans@gmail.com)

class: Workflow
cwlVersion: v1.0

doc: "jaffa workflow"

hints:
  - class: synData
    input: jaffa_reference
    entity: syn7885196

inputs:
  jaffa_reference:
    type: File

  REFERENCE_GENOME:
    type: File

  TUMOR_FASTQ_1:
    type: File

  TUMOR_FASTQ_2:
    type: File

outputs:

  OUTPUT:
    type: File
    outputSource: converttobedpe/fusionout

steps:
  tar:
    run: ../cwl/tar.cwl
    in:
      input: jaffa_reference
    out: [output]

  jaffa:
    run: ../cwl/jaffa.cwl
    in:
      threads: { default: 2 }
      pipeline: { default: "/opt/JAFFA/JAFFA_direct.groovy" }
      fastqInputFormat: { default: "%_*.fq.gz" }
      fastq1: TUMOR_FASTQ_1
      fastq2: TUMOR_FASTQ_2
      refBase: tar/output
      genomeFasta: REFERENCE_GENOME
      genome: { default: "GRCh37" }
      annotation: { default: "75" }
    out:
      [jaffa_fusions_list]

  converttobedpe:
    run: ../cwl/converter.cwl
    in:
      input: jaffa/jaffa_fusions_list
      output: { default: "output.bedpe" }
    out: [fusionout]

