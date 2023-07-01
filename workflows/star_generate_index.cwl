class: Workflow
cwlVersion: v1.0
id: star_generate_index
label: STAR_generate_index
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: reference
    type: File
    'sbg:x': -107
    'sbg:y': 27
  - id: threads
    type: int?
    'sbg:x': -12
    'sbg:y': -182
outputs:
  - id: star_index_directory
    outputSource:
      - generate_index/star_index_directory
    type: Directory
    'sbg:x': 534.7459106445312
    'sbg:y': 81.5
steps:
  - id: generate_index
    in:
      - id: reference
        source: reference
      - id: genome_chr_bin_n_bits
        source: calculate_star_index_generation_parameters/genomeChrBinNbits_file
      - id: threads
        source: threads
      - id: genome_sa_index_and_bases
        source: calculate_star_index_generation_parameters/genomeSAindexNbases_file
    out:
      - id: star_index_directory
    run: ../tools/star/2.7.10a/star_generate_index.cwl
    label: generate_index
    'sbg:x': 282.875
    'sbg:y': 10
  - id: calculate_star_index_generation_parameters
    in:
      - id: genome_reference_file
        source: reference
    out:
      - id: genomeSAindexNbases_file
      - id: genomeChrBinNbits_file
    run: ../tools/finder/2.0.0/calculate_star_index_generation_parameters.cwl
    label: calculate_STAR_indes_generation_parameters
    'sbg:x': 49.875
    'sbg:y': 135
requirements: []
