class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: generate_index
baseCommand:
  - STAR
inputs:
  - id: reference
    type: File
    inputBinding:
      position: 0
      prefix: '--genomeFastaFiles'
      shellQuote: false
  - id: genome_chr_bin_n_bits
    type: File
  - 'sbg:toolDefaultValue': '1'
    id: threads
    type: int?
    inputBinding:
      position: 0
      prefix: '--runThreadN'
      shellQuote: false
  - id: genome_sa_index_and_bases
    type: File
outputs:
  - id: star_index_directory
    type: Directory
    outputBinding:
      glob: '*star_index'
label: generate_index
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: '${return "--runMode genomeGenerate"}'
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: '${return "--genomeDir star_index"}'
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: >-
      ${return "--genomeSAindexNbases `cat "+
      inputs.genome_sa_index_and_bases.path + "`"}
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: >-
      ${return "--genomeChrBinNbits `cat
      "+inputs.genome_chr_bin_n_bits.path+"`"}
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    coresMin: 0
  - class: DockerRequirement
    dockerPull: 'sagnikbanerjee15/star:2.7.10a'
  - class: InlineJavascriptRequirement
stdout: star_index_generation.output
stderr: star_index_generation.error
