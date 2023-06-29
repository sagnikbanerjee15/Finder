# Proposed changes to FINDER

Here we list all the proposed changes that we intend to make to FINDER. The order of appearence does not reflect priority.

## Major architectural changes

- [ ] Change the directory structure to following:
	- [ ] src [will contain all the programs, including Code in C(if we decide to have some code)]
	- [ ] example [will contain some reference fasta files along with some fastq for testing installation]
	- [ ] dockerfiles [will contain dockerfiles for all the software involved. Each dockerfile needs to be in its own directory and within it will be a directory for the version. For example, tools/star/2.7.9a/Dockerfile]
	- [ ] tools [will contain CWL scripts for the tools. Each tool will be in its own directory and within it will be a directory for the version. For example, tools/star/2.7.9a]
	- [ ] workflows [will contain CWL scripts for workflows. All workflows must utilize tools and/or workflows from this repository alone.]
- [ ] Implement the entire workflow as a CWL pipeline - very easy to maintain and update software. cwltool, python and docker/singularity are the only software that needs to be installed on the parent system
- [ ] Provide option for a pip install and a conda install (have been a lot of requests from several users)
- [ ] BRAKER requires license files and makes it very difficult for installation. We will replace BRAKER with Augustus since our focus is more on RNA-Seq data rather than on predictions.
- [ ] FINDER2 will have a main program that will be named `finder`. This program will be written in python and will execute the workflows for various activities. This approach is favoured since a major part of the pipeline needs to execute on multiple samples and then aggregate the results into one and execute the rest of the pipeline
- [ ] Utmost care needs to be taken to ensure that there is minimal space consumption. For some datasets, especially those that have a very large number of samples, the space demand can skyrocket easily.

## The nitty gritties

- [ ] Decide where to have an ECR. Options are docker, ghcr
- [ ] Write up a program to enable users create docker/singularity images locally instead of pulling from any repo
- [ ] Implement the following Dockerfiles and tools
	- [ ] multiqc
	- [ ] repeatmasker (to perform soft masking to enhance mapping quality and detect genes that are present in low complexity regions) 
	- [ ] sratoolkit (to download data directly from NCBI SRA)
	- [ ] trimmomatic (for adapter trimming)
	- [ ] star (to perform alignment of short reads)
	- [ ] minimap2 (to perform alignment of long reads)
	- [ ] psiclass (to perform genome assembly)
	- [ ] codan (to perform coding region annotation)
	- [ ] Augustus (to predict genes)
	- [ ] changepoint in R (to split merged genes from RNA-Seq data)
	- [ ] exonerate (to align peptide sequences to reference genome)
	- [ ] NCBI-BLAST
	- [ ] regtools
	- [ ] samtools
	- [ ] bedtools
	- [ ] gffread
	- [ ] bedops
	- [ ] gffcompare
	- [ ] FINDER
	- [ ] Polyester (for simulating reads)
	- [ ] Other non coding RNA prediction tools 
	- [ ] DIAMOND
- [ ] Implement the following workflows
	- [ ] Short read alignment with star
	- [ ] Long read alignment with minimap2
	- [ ] Assembling short reads with psiclass
	- [ ] Execute codan for CDS prediction
	- [ ] Changepoint detection
	- [ ] Align peptide to reference genomes
- [ ] Decide on a proper format for accepting input from user. For this version, FINDER will be able to accept long reads data as well. Long reads transcriptomic data must be provided in fasta format.
- [ ] Collect data for a variety of species for testing FINDER2
- [ ] Decide on different metrics to judge quality of annotation for non-model organisms
- [ ] Decide on a proper method to accumulate outputs and errors from a run and send a summary to the developers. These could be the fields
	- sra_accession (SRA Accession number to download data directly from NCBI-SRA)
	- tissue_type (Type of tissue)
	- description (Free text - will not be used by FINDER, provided for user)
	- ended (PE or SE)
	- mrna_seq (Enter value 1 if data is mRNA-Seq, else provide 0 if data is small RNA-seq)
	- skip (by default FINDER will process all the samples, enter 1 if you wish to skip it but still retain the entry in the file. This option prevents the need to create additional metadata files)
	- location (Enter the location of the sample in the local disk)
	- same_species (Enter 1 if the RNA sample is from the same species. Else enter 0. Short read alignment will )
	- long_or_short
- [ ] Deploy final results through an HTML page. 
	- Integrate IGV and/or genome browsers for enhanced look and feel
	- Option to select RNA-Seq sample to view alignment
	- Provide option to select specific genes depending on several factors
- [ ] Output a single GTF file with more information about annotation source. Encode information in the "score" field of the GTF file only for transcripts. This is the 6th column.
	- Code represents a decimal number with each bit representing a special feature as listed below:

| Bit number   | Description |
|--------------|:-----|
| 0 | Predicted by AUGUSTUS & perfect overlap |
| 1 | Predicted by AUGUSTUS & partial overlap |
| 2 | Predicted by AUGUSTUS & no overlap |
| 3 | Peptide sequences provided |
| 4 | Found in peptide sequences |
| 5 | Has Coding sequence |
| 6 | Has both UTRs|
| 7 | Has only 5' UTR |
| 8 | Has only 3' UTR |
| 9 | Found in RNA-Seq data|
|10 | Has a single exon |
|11 | Constructed after changepoint split|
|12 | Represented in long reads|
|13 | Constructed after merging very closely placed transcript|
|14 | End exons adjusted after coverage|
|15 | Contains micro exons |
|16 | Is a small-RNA|
|17 | within 50 nucleotides of another transcript from a different gene on same strand|
|18 | within 50 nucleotides of another transcript from a different gene on opposite strand|
|19 | Contains introns larger than 10K|
|20 | Present in low complexity region|
|21 | Has less than 5 transcripts|
|22 | Has >5 but <=10 transcripts|
|23 | Has >10 but <=20 transcripts|
|24 | Has >20 transcripts|
|25-| Is present in xyz tissue type where xyz is a tissue type


- [ ] Output a peptide file 
- [ ] Deploy a mongodb database for storing genes generated at each phase
- [ ] Short read alignment - perform 3 rounds
- [ ] Configure sra-toolkit to download fasta files instead of fastq to save space
- [ ] Allow processing of small RNA-Seq data to annotate a variety of non coding RNAs
- [ ] Write a program to select different types of genes from the GTF file produced by FINDER
- [ ] Option to recreate gene annotations over existing gene annotations using new data (either long or short)
- [ ] Select the model and non-model organisms on which FINDER2 will be tested. Also select a wide range of RNA-Seq samples and tissue types. BRAKER3 is coming out so we can test with that as well. 
