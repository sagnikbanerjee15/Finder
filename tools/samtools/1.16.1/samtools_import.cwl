class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: samtools_import
baseCommand:
  - samtools import
inputs:
  - id: mate1
    type: File
  - id: mate2
    type: File?
outputs:
  - id: bam
    type: File?
    outputBinding:
      glob: '*bam'
label: samtools_import
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          if(inputs.mate2)
          {
              return "-1 "+inputs.mate1.path + " -2 "+inputs.mate2.path
          }
          else
          {
              return "-0 "+inputs.mate1.path
          }
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "-o "+" unaligned.bam"
      }
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'sagnikbanerjee15/samtools:1.16.1'
