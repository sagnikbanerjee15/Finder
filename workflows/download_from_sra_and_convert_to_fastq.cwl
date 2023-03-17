class: Workflow
cwlVersion: v1.0
id: download_from_sra_and_convert_to_fastq
label: download_from_sra_and_convert_to_fastq
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: threads
    type: int?
    'sbg:x': -616
    'sbg:y': -74
  - id: SRA_accession
    type: string
    'sbg:x': -626
    'sbg:y': 182
outputs:
  - id: pair2_fasta
    outputSource:
      - sratools_fasterq_dump/pair2_fasta
    type: File?
    'sbg:x': -218.109375
    'sbg:y': -123.5
  - id: pair1_fasta
    outputSource:
      - sratools_fasterq_dump/pair1_fasta
    type: File?
    'sbg:x': -120.109375
    'sbg:y': 107.5
  - id: fastas
    outputSource:
      - sratools_fasterq_dump/fastas
    type: 'File[]'
    'sbg:x': -238.109375
    'sbg:y': 327.5
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
    run: ../tools/sratools/3.0.0/fasterq_dump.cwl
    label: sratools_fasterq_dump
    'sbg:x': -350
    'sbg:y': 46
requirements: []
