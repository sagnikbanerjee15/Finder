class: Workflow
cwlVersion: v1.0
id: align_with_star_multiple_rounds
label: align_with_star_multiple_rounds
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: star_index
    type: Directory
    'sbg:x': 0
    'sbg:y': 137.609375
  - id: threads
    type: int?
    'sbg:x': 0
    'sbg:y': 30.390625
  - id: raw_input_files_pair2
    type: File?
    'sbg:x': 0
    'sbg:y': 244.828125
  - id: raw_input_files_pair1
    type: File
    'sbg:x': 0
    'sbg:y': 352.046875
  - id: max_multimaps_round_1
    type: int?
    'sbg:x': 0
    'sbg:y': 888.140625
  - id: min_overhang_length_for_SJ_in_DB_round_1
    type: int?
    'sbg:x': 0
    'sbg:y': 780.921875
  - id: min_overhang_length_for_SJ_not_in_DB_round_1
    type: int?
    'sbg:x': 0
    'sbg:y': 673.703125
  - id: max_multimaps_round_2
    type: int?
    'sbg:x': 456.6875
    'sbg:y': 727.3125
  - id: min_overhang_length_for_SJ_not_in_DB_round_2
    type: int?
    'sbg:x': 456.6875
    'sbg:y': 512.875
  - id: min_overhang_length_for_SJ_in_DB_round_2
    type: int?
    'sbg:x': 456.6875
    'sbg:y': 620.09375
  - id: min_overhang_length_for_SJ_not_in_DB_round_3
    type: int?
    'sbg:x': 1630.5103759765625
    'sbg:y': 1000
  - id: min_overhang_length_for_SJ_in_DB_round_3
    type: int?
    'sbg:x': 1635.6083984375
    'sbg:y': 1118.19580078125
  - id: max_multimaps_round_3
    type: int?
    'sbg:x': 1639.524169921875
    'sbg:y': 1242.6275634765625
  - id: outFilterMismatchNoverLmax_round1
    type: float?
    'sbg:x': 0
    'sbg:y': 459.265625
  - id: outFilterMismatchNoverLmax_round2
    type: float?
    'sbg:x': 456.6875
    'sbg:y': 298.4375
  - id: outFilterMismatchNoverLmax_round3
    type: float?
    'sbg:x': 1620.9229736328125
    'sbg:y': 765.80419921875
  - id: outFilterMatchNminOverLread_round1
    type: float?
    'sbg:x': 0
    'sbg:y': 566.484375
  - id: outFilterMatchNminOverLread_round2
    type: float?
    'sbg:x': 456.6875
    'sbg:y': 405.65625
  - id: outFilterMatchNminOverLread_round3
    type: float?
    'sbg:x': 1627.0208740234375
    'sbg:y': 878.19580078125
outputs:
  - id: log_round_1
    outputSource:
      - star_align_round_1/log
    type: File
    'sbg:x': 898.2645263671875
    'sbg:y': -50.6320915222168
  - id: log_round_2
    outputSource:
      - star_align_round_2/log
    type: File
    'sbg:x': 1706.5704345703125
    'sbg:y': -30.233205795288086
  - id: log_round_3
    outputSource:
      - star_align_round_3/log
    type: File
    'sbg:x': 2349.228759765625
    'sbg:y': 501.32647705078125
  - id: merged_bam
    outputSource:
      - samtools_merge/output_bam
    type: File?
    'sbg:x': 2664.881591796875
    'sbg:y': 1042.9794921875
  - id: star_align_round1_output
    outputSource:
      - star_align_round_1/star_align_output
    type: File
    'sbg:x': 896.3431396484375
    'sbg:y': -290.0250549316406
  - id: star_align_round1_error
    outputSource:
      - star_align_round_1/star_align_error
    type: File
    'sbg:x': 899.4476928710938
    'sbg:y': -167.18409729003906
  - id: star_align_round3_output
    outputSource:
      - star_align_round_3/star_align_output
    type: File
    'sbg:x': 2337.6318359375
    'sbg:y': 259.6318054199219
  - id: star_align_round3_error
    outputSource:
      - star_align_round_3/star_align_error
    type: File
    'sbg:x': 2342.3681640625
    'sbg:y': 391.3431396484375
  - id: star_align_round2_output
    outputSource:
      - star_align_round_2/star_align_output
    type: File
    'sbg:x': 1696.6318359375
    'sbg:y': -272.0250549316406
  - id: star_align_round2_error
    outputSource:
      - star_align_round_2/star_align_error
    type: File
    'sbg:x': 1703.41455078125
    'sbg:y': -154.26419067382812
steps:
  - id: star_align_round_1
    in:
      - id: threads
        source: threads
      - id: genome_directory
        source: star_index
      - id: raw_input_files_pair1
        source: raw_input_files_pair1
      - id: raw_input_files_pair2
        source: raw_input_files_pair2
      - id: min_intron_length
        default: 20
      - id: max_intron_length
        default: 10000
      - id: RG_id
        default: 1
      - id: sam_attributes
        default: All
      - id: max_multimaps
        source: max_multimaps_round_1
      - id: max_memory
        default: 50
      - id: allow_soft_clipping
        default: true
      - id: min_overhang_length_for_SJ_not_in_DB
        source: min_overhang_length_for_SJ_not_in_DB_round_1
      - id: min_overhang_length_for_SJ_in_DB
        source: min_overhang_length_for_SJ_in_DB_round_1
      - id: outFilterMismatchNoverLmax
        source: outFilterMismatchNoverLmax_round1
      - id: outFilterMatchNminOverLread
        source: outFilterMatchNminOverLread_round1
    out:
      - id: alignment_file
      - id: log
      - id: splice_junctions
      - id: unmapped_1
      - id: unmapped_2
      - id: star_align_output
      - id: star_align_error
    run: ../tools/star/2.7.10a/star_align.cwl
    label: star_align_round_1
    'sbg:x': 456.6875
    'sbg:y': 135.21875
  - id: star_align_round_2
    in:
      - id: genome_directory
        source: star_index
      - id: raw_input_files_pair1
        source: star_align_round_1/unmapped_1
      - id: raw_input_files_pair2
        source: star_align_round_1/unmapped_2
      - id: min_intron_length
        default: 20
      - id: max_intron_length
        default: 10000
      - id: RG_id
        default: 2
      - id: sam_attributes
        default: All
      - id: max_multimaps
        source: max_multimaps_round_2
      - id: max_memory
        default: 50
      - id: allow_soft_clipping
        default: true
      - id: min_overhang_length_for_SJ_not_in_DB
        source: min_overhang_length_for_SJ_not_in_DB_round_2
      - id: min_overhang_length_for_SJ_in_DB
        source: min_overhang_length_for_SJ_in_DB_round_2
      - id: splice_junction_db_file
        source: star_align_round_1/splice_junctions
      - id: outFilterMismatchNoverLmax
        source: outFilterMismatchNoverLmax_round2
      - id: outFilterMatchNminOverLread
        source: outFilterMatchNminOverLread_round2
    out:
      - id: alignment_file
      - id: log
      - id: splice_junctions
      - id: unmapped_1
      - id: unmapped_2
      - id: star_align_output
      - id: star_align_error
    run: ../tools/star/2.7.10a/star_align.cwl
    label: star_align_round_2
    'sbg:x': 1381.6171875
    'sbg:y': 107.10894775390625
  - id: star_align_round_3
    in:
      - id: threads
        source: threads
      - id: genome_directory
        source: star_index
      - id: raw_input_files_pair1
        source: star_align_round_2/unmapped_1
      - id: raw_input_files_pair2
        source: star_align_round_2/unmapped_2
      - id: min_intron_length
        default: 10000
      - id: max_intron_length
        default: 100000
      - id: RG_id
        default: 3
      - id: sam_attributes
        default: All
      - id: max_multimaps
        source: max_multimaps_round_3
      - id: max_memory
        default: 50
      - id: allow_soft_clipping
        default: true
      - id: min_overhang_length_for_SJ_not_in_DB
        source: min_overhang_length_for_SJ_not_in_DB_round_3
      - id: min_overhang_length_for_SJ_in_DB
        source: min_overhang_length_for_SJ_in_DB_round_3
      - id: outFilterMismatchNoverLmax
        source: outFilterMismatchNoverLmax_round3
      - id: outFilterMatchNminOverLread
        source: outFilterMatchNminOverLread_round3
      - id: seedPerWindowNmax
        default: 10
    out:
      - id: alignment_file
      - id: log
      - id: splice_junctions
      - id: unmapped_1
      - id: unmapped_2
      - id: star_align_output
      - id: star_align_error
    run: ../tools/star/2.7.10a/star_align.cwl
    label: star_align_round_3
    'sbg:x': 2073.285888671875
    'sbg:y': 455.20745849609375
  - id: samtools_merge
    in:
      - id: input_alignment
        source:
          - star_align_round_1/alignment_file
          - star_align_round_2/alignment_file
          - star_align_round_3/alignment_file
      - id: output_format
        default: BAM
      - id: threads
        source: threads
    out:
      - id: output_bam
      - id: output_sam
      - id: output_cram
      - id: stdout
      - id: stderr
    run: ../tools/samtools/1.16.1/samtools_merge.cwl
    label: samtools merge
    'sbg:x': 2364.959228515625
    'sbg:y': 831.4197998046875
requirements:
  - class: MultipleInputFeatureRequirement
'sbg:toolAuthor': Sagnik Banerjee
