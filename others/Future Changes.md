# Changelog

#### High Priority (Will be addressed in the next release)

- Create a docker image
- Removed BRAKER2 implementation
- Removed GeneMarkS/T for predicting coding sequence
- Improved gene annotation by adding more information about the source of evidence RNA-Seq, protein and/or prediction. Reduced the number of GTF output files.
- Checks for low ulimit
- Remove the option to include path locations for individual files and modify the FINDER to accept a single path. The entire directory tree will be searched. Repeated occurance will generate error.
- Added a functionality to predict genes using AUGUSTUS (https://math-inf.uni-greifswald.de/storages/uni-greifswald/fakultaet/mnf/mathinf/stanke/augustus_wrp.pdf)
- Added a new software `CodAn` to predict coding sequence
- Added an advanced option to include a fraction to denote the maximum number of repeat elements for a transcript/gene to be considered low confidence
- Update all tools to latest version 
- Add an option to decide the overlap of braker transcripts with proteins


#### Lower priority

- Improve logging for each submodule such that all errors can be easily spotted
- Update the program to download fastq files. Now it will download gzipped files to reduce space demands
- Change the implementation of Change Point Detection (CPD). Rewritte the code in C to speed up execution
- Change the method of alignment. Combined regions of genome with low expression and assemble them separately
- Change BLAST to DIAMOND to speed up the pipeline
- Change the method to generate alignments to micro-exons. Instead of mapping all unmapped reads using OLEGO, we now target reads that were soft-clipped around the presumed micro-exon
- Adde a program `xyz` to select/remove genes/transcripts. Also it will report the genes that overlap with one another (protein-protein, protein-ncRNA, ncRNA-ncRNA).
- Adde a functionality to include predictions from multiple sources
- Incorporate ABRIDGE to compress intermediate alignment files
- Add functionality to process long read data from PacBio and/or Nanopore
- Add RepeatMasker for genomes
- Add the functionality to look for fastq, fasta, fasta.gz, fastq.gz, fq, fq.gz, fa, fa.gz
- Add the option to adding memory constraint
- Add an option to include an existing annotation (as a GTF file). FINDER will only generate and report novel gene models
- Add options to process small RNA-Seq data and annotate non-coding transcripts
- Write a program to collect all logs, output and error files and generate a zipped file that can be emailed directly to the developers.
- Add the functionality to process CAGE-Seq, RAMPAGE-Seq reads
- Output gene fasta and protein fasta files
- Add functionality to include gene models from closely related species
- Add functionality to annotate proteins with GO terms or with domain information
- Output a summary of the different gene structures like number of genes from prediction, no. of genes from expression data only etc.
- Include output in GFF3 with all information like UTR, start and stop codons
- Support for processing strand-specific RNA-Seq data
- Include BUSCO to quality check predicted transcripts
- Remove concept of low or high confidence genes/transcripts. Put out all the levels of information that was used to construct it.
- Integrate genomebrowser and allow to view big wig files along with annotations [major change]
- Merge contigs and split them to speed up processing
- Parallelize over several nodes - useful for HPCs
- Group transcripts into correct set of genes. 
- Modify SAM files to contain only mappings - no sequence and no quality scores