class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: samtools_fqidx
baseCommand:
  - samtools
  - fqidx
inputs:
  - id: input_fastq
    type: File
    inputBinding:
      position: 100
      shellQuote: false
outputs:
  - id: output
    type: File?
    outputBinding:
      glob: '*fai'
  - id: stdout
    type: File?
    outputBinding:
      glob: '*output'
  - id: stderr
    type: File?
    outputBinding:
      glob: '*error'
label: samtools_fqidx
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'sagnikbanerjee15/samtools:1.16.1'
stdout: samtools_fqids.output
stderr: samtools_fqids.error
