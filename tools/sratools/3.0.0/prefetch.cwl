class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: sratools_prefetch
baseCommand:
  - prefetch
inputs:
  - id: SRA_accession
    type: string
    inputBinding:
      position: 100
      shellQuote: false
outputs:
  - id: sra_file
    type: File
    outputBinding:
      glob: '*sra'
label: sratools_prefetch
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--type sra"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--output-directory ."
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--quiet"
      }
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ncbi/sra-tools:3.0.1'
