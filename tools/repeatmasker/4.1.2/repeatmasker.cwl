class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: repeatmasker
baseCommand:
  - RepeatMasker
inputs:
  - id: sequence_fasta_file
    type: File
    inputBinding:
      position: 100
      shellQuote: false
  - id: threads
    type: int
    inputBinding:
      position: 0
      prefix: '-pa'
      shellQuote: false
  - id: species
    type:
      type: enum
      symbols:
        - human
        - mouse
        - rattus
        - arabidopsis
        - mammal
        - carnivore
        - rodentia
        - rat
        - cow
        - pig
        - cat
        - dog
        - chicken
        - fugu
        - danio
        - anopheles
        - worm
        - diatoaea
        - artiodactyl
        - rice
        - wheat
        - maize
        - algae
      name: species
    inputBinding:
      position: 0
      prefix: '--species'
      shellQuote: false
outputs:
  - id: masked_fasta
    type: File
    outputBinding:
      glob: '*_repeatmasker/*masked'
label: repeatmasker
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "-qq"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--norna"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "-dir "+inputs.sequence_fasta_file.nameroot+"_repeatmasker"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "-small -xsmall -lcambig"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "-html"
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'sagnikbanerjee15/repeatmasker:4.1.2'
  - class: InlineJavascriptRequirement
