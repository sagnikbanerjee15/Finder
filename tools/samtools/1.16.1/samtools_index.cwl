class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: samtools_index
baseCommand: []
inputs:
  - id: bamfilename
    type: File
    doc: Bam file name
  - 'sbg:toolDefaultValue': '1'
    id: threads
    type: int?
  - id: cai_index
    type: boolean?
outputs:
  - id: bai
    type: File?
    outputBinding:
      glob: '*.bai'
  - id: cai
    type: File?
    outputBinding:
      glob: '*cai'
  - id: stdout
    type: File?
    outputBinding:
      glob: '*.output'
  - id: stderr
    type: File?
    outputBinding:
      glob: '*.error'
label: samtools_index
arguments:
  - position: 6
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "-@ " + inputs.threads
      }
  - position: 8
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return " " + inputs.bamfilename.basename
          
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "cp "+inputs.bamfilename.path+" . && samtools index "
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'sagnikbanerjee15/samtools:1.16.1'
  - class: InlineJavascriptRequirement
stdout: samtools_index.output
stderr: samtools_index.error
