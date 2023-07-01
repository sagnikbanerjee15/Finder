class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: calculate_star_index_generation_parameters
baseCommand:
  - python
  - calculate_STAR_index_generation_parameters
inputs:
  - id: genome_reference_file
    type: File
    inputBinding:
      position: 0
      prefix: '--reference'
      shellQuote: false
outputs:
  - id: genomeSAindexNbases_file
    type: File
    outputBinding:
      glob: '*genomeSAindexNbases'
  - id: genomeChrBinNbits_file
    type: File
    outputBinding:
      glob: '*genomeChrBinNbits'
label: calculate_STAR_indes_generation_parameters
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--genomeSAindexNbases_filename "+inputs.genome_reference_file.nameroot+"_genomeSAindexNbases"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--genomeChrBinNbits_filename "+inputs.genome_reference_file.nameroot+"_genomeChrBinNbits"
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'sagnikbanerjee15/star:2.7.10a'
  - class: InlineJavascriptRequirement
