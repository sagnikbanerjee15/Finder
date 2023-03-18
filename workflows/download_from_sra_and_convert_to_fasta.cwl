class: Workflow
cwlVersion: v1.0
id: download_from_sra_and_convert_to_fastq
label: download_from_sra_and_convert_to_fastq
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: threads
    type: int?
    'sbg:x': 0
    'sbg:y': 160.5
  - id: SRA_accession
    type: string
    'sbg:x': 0
    'sbg:y': 267.5
outputs:
  - id: pair2_fasta
    outputSource:
      - sratools_fasterq_dump/pair2_fasta
    type: File?
    'sbg:x': 513.29345703125
    'sbg:y': 0
  - id: pair1_fasta
    outputSource:
      - sratools_fasterq_dump/pair1_fasta
    type: File?
    'sbg:x': 513.29345703125
    'sbg:y': 107
  - id: fastas
    outputSource:
      - sratools_fasterq_dump/fastas
    type: 'File[]'
    'sbg:x': 513.29345703125
    'sbg:y': 428
  - id: fasterq_dump_output
    outputSource:
      - sratools_fasterq_dump/fasterq_dump_output
    type: File
    'sbg:x': 513.29345703125
    'sbg:y': 214
  - id: fasterq_dump_error
    outputSource:
      - sratools_fasterq_dump/fasterq_dump_error
    type: File
    'sbg:x': 513.29345703125
    'sbg:y': 321
steps:
  - id: sratools_fasterq_dump
    in:
      - id: SRA_accession
        source: SRA_accession
      - id: threads
        source: threads
    out:
      - id: fastas
      - id: pair1_fasta
      - id: pair2_fasta
      - id: fasterq_dump_output
      - id: fasterq_dump_error
    run: ../tools/sratools/3.0.0/fasterq_dump.cwl
    label: sratools_fasterq_dump
    'sbg:x': 175.453125
    'sbg:y': 186
requirements: []
