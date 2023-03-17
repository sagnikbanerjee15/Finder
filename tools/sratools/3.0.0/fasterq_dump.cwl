class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: sratools_fasterq_dump
baseCommand:
  - fasterq-dump
inputs:
  - id: SRA_accession
    type: string
    inputBinding:
      position: 100
      shellQuote: false
  - 'sbg:toolDefaultValue': '1'
    id: threads
    type: int?
    inputBinding:
      position: 0
      prefix: '--threads'
      shellQuote: false
outputs:
  - id: fastas
    type: 'File[]'
    outputBinding:
      glob: '*fasta'
  - id: pair1_fasta
    type: File?
    outputBinding:
      glob: '*_1.fasta'
  - id: pair2_fasta
    type: File?
    outputBinding:
      glob: '*_2.fasta'
label: sratools_fasterq_dump
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--progress"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--outdir ."
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--split-files"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--seq-defline '>$ac_$sn_$sg_$si_$rl/$ri'"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--fasta"
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ncbi/sra-tools:3.0.1'
  - class: InlineJavascriptRequirement
