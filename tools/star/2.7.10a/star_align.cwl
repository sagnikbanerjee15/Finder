class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: star_align
baseCommand: []
inputs:
  - 'sbg:toolDefaultValue': '1'
    id: threads
    type: int?
    inputBinding:
      position: 0
      prefix: '--runThreadN'
      shellQuote: false
  - id: genome_directory
    type: Directory
  - id: raw_input_files_pair1
    type: File
  - id: raw_input_files_pair2
    type: File?
  - 'sbg:toolDefaultValue': '20'
    id: min_intron_length
    type: int?
    inputBinding:
      position: 0
      prefix: '--alignIntronMin'
      shellQuote: false
  - 'sbg:toolDefaultValue': '10000'
    id: max_intron_length
    type: int?
    inputBinding:
      position: 0
      prefix: '--alignIntronMax'
      shellQuote: false
  - id: RG_id
    type: int?
  - id: sam_attributes
    type: string?
    inputBinding:
      position: 0
      prefix: '--outSAMattributes'
      shellQuote: false
  - 'sbg:toolDefaultValue': '500'
    id: max_multimaps
    type: int?
    inputBinding:
      position: 0
      prefix: '--outFilterMultimapNmax'
      shellQuote: false
  - 'sbg:toolDefaultValue': '1'
    id: max_memory
    type: int?
  - id: allow_soft_clipping
    type: boolean
  - 'sbg:toolDefaultValue': '8'
    id: min_overhang_length_for_SJ_not_in_DB
    type: int?
    inputBinding:
      position: 0
      prefix: '--alignSJoverhangMin'
      shellQuote: false
  - 'sbg:toolDefaultValue': '8'
    id: min_overhang_length_for_SJ_in_DB
    type: int?
    inputBinding:
      position: 0
      prefix: '--alignSJDBoverhangMin'
      shellQuote: false
  - id: splice_junction_db_file
    type: File?
    inputBinding:
      position: 0
      prefix: '--sjdbFileChrStartEnd'
      shellQuote: false
  - 'sbg:toolDefaultValue': '0.3'
    id: outFilterMismatchNoverLmax
    type: float?
    inputBinding:
      position: 0
      prefix: '--outFilterMismatchNoverLmax'
      shellQuote: false
  - 'sbg:toolDefaultValue': '0.66'
    id: outFilterMatchNminOverLread
    type: float?
    inputBinding:
      position: 0
      prefix: '--outFilterMatchNminOverLread'
      shellQuote: false
outputs:
  - id: alignment_file
    type: File
    outputBinding:
      glob: '*bam'
  - id: log
    type: File
    outputBinding:
      glob: '*log'
  - id: splice_junctions
    type: File
    outputBinding:
      glob: '*.tab'
  - id: unmapped_1
    type: File?
    outputBinding:
      glob: '*unmapped*1*fastq'
  - id: unmapped_2
    type: File?
    outputBinding:
      glob: '*unmapped*2*fastq'
label: star_align
arguments:
  - position: 1
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--outSAMtype BAM SortedByCoordinate"
      }
  - position: 2
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--outReadsUnmapped Fastx"
      }
  - position: 3
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          if(inputs.RG_id)
          {
              return "--outSAMattrRGline ID:"+inputs.RG_id
          }
      }
  - position: 4
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--outFileNamePrefix _"
      }
  - position: 5
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          if(inputs.allow_soft_clipping==false)
          {
              return "--alignEndsType EndToEnd"
          }
          else
          {
              return "--alignEndsType Local"
          }
      }
  - position: 6
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          if(inputs.raw_input_files_pair2 == undefined)
          {
              return "--readFilesIn " + inputs.raw_input_files_pair1.path
          }
          else
          {
              return "--readFilesIn " + inputs.raw_input_files_pair1.path + " "+ inputs.raw_input_files_pair2.path
          }
      }
  - position: 7
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          return "--limitBAMsortRAM "+inputs.max_memory * 1024 * 1024 * 1024
      }
  - position: 100
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          var cmd1 = "mv _Aligned.sortedByCoord.out.bam " + inputs.raw_input_files_pair1.nameroot + "_possorted.bam"
          
          var cmd2 = "mv _Unmapped.out.mate1 " + inputs.raw_input_files_pair1.nameroot + "_unmapped_1.fastq"
          
          
          
          var cmd4 = "mv _Log.final.out " + inputs.raw_input_files_pair1.nameroot + ".log"
          
          var cmd5 = "mv _SJ.out.tab " + inputs.raw_input_files_pair1.nameroot + ".tab"
          
          if(inputs.raw_input_files_pair2 == null)
          {
              return "&& " + cmd1 + " && " + cmd2 + " && " + cmd4 + " && " + cmd5
          }
          else
          {
              var cmd3 = "mv _Unmapped.out.mate2 " + inputs.raw_input_files_pair2.nameroot + "_unmapped_2.fastq"
              return "&& " + cmd1 + " && " + cmd2 + " && " + cmd3 + " && " + cmd4 + " && " + cmd5
          }
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          if(inputs.splice_junction_db_file == undefined)
              return "STAR --genomeDir "+inputs.genome_directory.path
          else
              return "mkdir star_index && cp -r "+inputs.genome_directory.path+"/* star_index/ && STAR --genomeDir star_index"
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'sagnikbanerjee15/star:2.7.10a'
  - class: InlineJavascriptRequirement
