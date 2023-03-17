class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: samtools_faidx
baseCommand:
  - ln
  - '-s'
inputs:
  - id: fa
    type: File
    inputBinding:
      position: 1
    doc: FastA file for reference genome
outputs:
  - id: fai
    type: File
    outputBinding:
      glob: '*.fai'
  - id: stdout
    type: File?
    outputBinding:
      glob: '*.output'
  - id: stderr
    type: File?
    outputBinding:
      glob: '*.error'
label: samtools_faidx
arguments:
  - position: 2
    prefix: ''
    shellQuote: false
    valueFrom: $(inputs.fa.basename)
  - position: 3
    prefix: ''
    shellQuote: false
    valueFrom: '&&'
  - position: 4
    prefix: ''
    shellQuote: false
    valueFrom: samtools
  - position: 5
    prefix: ''
    shellQuote: false
    valueFrom: faidx
  - position: 6
    prefix: ''
    shellQuote: false
    valueFrom: $(inputs.fa.basename)
  - position: 7
    prefix: ''
    shellQuote: false
    valueFrom: '&&'
  - position: 8
    prefix: ''
    shellQuote: false
    valueFrom: rm
  - position: 9
    prefix: ''
    shellQuote: false
    valueFrom: $(inputs.fa.basename)
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'sagnikbanerjee15/samtools:1.16.1'
stdout: samtools_faidx.output
stderr: samtools_faidx.error
