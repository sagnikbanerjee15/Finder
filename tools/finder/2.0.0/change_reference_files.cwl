class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: change_reference_files
baseCommand:
  - python
  - change_reference_ids
inputs:
  - id: reference
    type: 'File[]'
    inputBinding:
      position: 0
      prefix: '--reference'
      shellQuote: false
outputs:
  - id: modified_reference
    type: File
    outputBinding:
      glob: '*modified_reference.fasta'
label: change_reference_files
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--mapping change_reference_ids.log"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--modified_reference modified_reference.fasta"
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ghcr.io/tools_and_pipelines/finder:2.0.0'
  - class: InlineJavascriptRequirement
