# Welcome to FINDER

FINDER is a gene annotator pipeline which automates the process of downloading short reads, aligning them and using the assembled  transcripts to generate gene annotations. Additionally it used protein sequences and reports gene predictions by BRAKER2. It is a fast, scalable, platform independent software that generatess gene annotations in GTF format. FINDER accepts inputs through command line interface. It finds several novel genes/transcripts and also reports the tissue/conditions they were found to be in. If you use FINDER for your research please cite <>

## Installation

FINDER requires a number of softwares which needs to be installed. This might cause version conflicts with softwares that are already installed in your system. Hence, the developers have decided to enforce the use of FINDER within a conda environment. 

### Downloading and installing Anaconda (skip if you already have Anaconda installed)

Execute the following commands to download and install Anaconda

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
bash Anaconda3-2020.02-Linux-x86_64.sh # do not install VS. You may replace it with the latest version of Anaconda from https://www.anaconda.com/distribution/. Also, conda will default to the home directory, but make sure you choose a directory that has sufficient disk space to install all the software packages. 
```

## Downloading FINDER and creating a conda environment

```bash
git clone https://github.com/sagnikbanerjee15/finder.git
cd finder
conda env create -f environment.yml # This will create an environment named finder_conda_env
conda activate finder
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc # Add this path permanently to the bashrc file
export PATH=$PATH:$(pwd)
cd dep
```

FINDER runs BRAKER2 which depends on GeneMark-ES. FINDER also needs GeneMarkS/T to predict coding regions of genes. Both GeneMark-ES and GeneMarkS/T are hosted at the University of Georgia website. The license prohibits the redistribution of their software, which is why it could not be included in this package. Hence, users have to manually download these 2 softwares and place them under the src sub-directory. Please follow the instructions below to download the softwares and the key:

1. Open a browser of your choice
2. Go to [this](http://topaz.gatech.edu/GeneMark/license_download.cgi) website
3. Select the option **GeneMark-ES/ET/EP ver 4.62_lic** (2<sup>nd</sup> from top) and **LINUX 64**
4. Enter your name, institution, country and email-id and click on the button that says ***I agree to the terms of this license agreement****
5. Right click on the link that says *Please download program **here*** and select **Copy Link Address**
6. Then type in `wget ` and paste the path you just copied
7. This command will download the file **gmes_linux_64.tar.gz** in the `src` sub-directory
8. Go to step 2
9. Select the option **GeneMarkS-T** (last option) and **LINUX 64**
10. Repeat steps 4-6
11. This command will download the file **gmst_linux_64.tar.gz** in the `src` sub-directory
12. Now, right click on the link that says **64_bit** and select **Copy Link Address**
13. Then type in `wget ` and paste the path you just copied
14. This command will download the file **gm_key_64.tar.gz** in the `src` sub-directory. This license key will serve for both the programs. Please note that this key will expire after one year from the date of download.
15. Execute the following commands:

```bash
./install.sh
```



## Executing FINDER with Sample data

Please follow the following the instructions to generate gene annotations using *Arabidopsis thaliana*. A CSV file template has been provided with the release in `example/Arabidopsis_thaliana_metadata.csv`. Keep all the headers intact and replace the data with your samples of choice. Also note, that FINDER can work with both data downloaded from NCBI and also with data on local directories. Below is a detailed description of the each column of the metadata file.

| Column Name      | Column Description                                           | Mandatory |
| :--------------- | :----------------------------------------------------------- | :-------- |
| BioProject       | Name of the bioproject that the data belongs to. If you are using locally saved data then please enter a dummy project name. Please note that FINDER will **NOT** be able to process empty fields of Bioproject. | **YES**   |
| SRA Accession    | Enter the SRA Accession number of the sample that you expect FINDER to use for generating the gene annotations. Note that FINDER will use this ID to download the read samples from NCBI-SRA. In case you wish to use data which is not currently uploaded to NCBI, then you should enter the name of the local file. Do not enter any file extension in this field. For example, in the template, we have used 6 local RNA-Seq samples for annotation. On line 19, we fill this cell with Atx_tm_A_rep_1 since it is the sample name. FINDER expects a single file for single-ended data and two files for paired-ended data. Single ended files must be named as <SRA_Accession>.fastq and the paired ended files must be named as <SRA_Accession>\_1.fastq and <SRA_Accession>\_2.fastq | **YES**   |
| Tissues          | Mention the tissue type or condition from which the sample has been collected. FINDER will report the tissues that are associated with a particular transcript. This can be used to find gene models that are expressed in a specific tissue and/or condition | **YES**   |
| Description      | A brief description of the data. This field is not mandatory and is not used by FINDER. It is upto the user to enter whatever metadata is deemed important. | NO        |
| Date             | Enter the date of producing the RNA-Seq sample. This field is not mandatory and is not used by FINDER. | NO        |
| Read Length (bp) | Enter the length of the reads. This field is not mandatory and is not used by FINDER. | NO        |
| Ended            | Enter either PE or SE for Paired ended reads or single neded reads. No other value should be entered. | **YES**   |
| RNA-Seq          | Enter 1 for all the rows. This field is included for future extensions. | **YES**   |
| process          | Enter 1 if you wish to process the sample. If a value of 0 is present, then FINDER will ignore the sample | **YES**   |
| Location         | Enter the location of the directory. For samples to be downloaded from NCBI, this field should be left empty. If the location of a directory is provided here then FINDER will assume that the sample is present in it. FINDER will generate an error if the sample is not found in this directory. It is *not* necessary to have all the samples in the same directory. To illustrate this, the template has 6 samples located in diferent folder locations. If you want to use pre-downloaded datat which is also available in NCBI, just provide the directory where the samples can be found. FINDER is configured to skip downloading from NCBI, if a local directory is specified. | **YES**   |

To optimize disk space usage FINDER will process read samples from each bioproject at a time. Once the data is downloaded and reads are mapped, FINDER will remove all those data (if `-no-cleanup` is not specificied) to save disk space. But samples that were locally present will not be removed.



Restart computations

Progress.log



Restart operations after failure

Restarting previous runs with more RNA-Seq samples

### Preparing the genome index

FINDER uses STAR and OLego for aligning and PsiCLASS for assembling. Users have the choice of generating the genome index and providing it to FINDER or allowing FINDER to generate the index. For some organisms with large genomes, STAR might require more memory. Hence users might be forced to transfer the STAR genome index generated in a computer that permits the usage of a large enough memory for index generation.  Download the *Arabidopsis thaliana* genome using the following command:

```bash
cd example
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
```

Unzip the genome using this command

```bash
gunzip Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
```

This command will generate the file `Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa` which contains the entire genome of *Arabidopsis thaliana*

STAR genome index can be created by the following command

```bash
CPU=30 #Enter the number of CPUs that you are permitted to use
mkdir star_index_without_transcriptome

STAR --runMode genomeGenerate --runThreadN $CPU --genomeDir star_index_without_transcriptome --genomeSAindexNbases 12 --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa
```

Similarly olego index can be generated by the following command

```bash
../dep/olego/olegoindex -p olego_index Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa
```

### Running FINDER

Help menu for FINDER can be launched by the following command:

```
FINDER -h

usage: FINDER [-h] --metadatafile METADATAFILE --output_directory
              OUTPUT_DIRECTORY --genome GENOME [--cpu CPU]
              [--genome_dir_star GENOME_DIR_STAR]
              [--genome_dir_olego GENOME_DIR_OLEGO] [--verbose VERBOSE]
              [--protein PROTEIN] [--no_cleanup] [--preserve_raw_input_data]
              [--checkpoint CHECKPOINT]

Generates gene annotation from RNA-Seq data

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  --metadatafile METADATAFILE, -mf METADATAFILE
                        Please enter the name of the metadata file. Enter 0 in the last column of those samples which you wish to skip processing. The columns should represent the following in order --> BioProject, SRA Accession, Tissues, Description, Date, Read Length, Ended (PE or SE), RNA-Seq, process, Location. If the sample is skipped it will not be downloaded. Leave the directory path blank if you are downloading the samples. In the end of the run the program will output a csv file with the directory path filled out. Please check the provided csv file for more information on how to configure the metadata file. 
  --output_directory OUTPUT_DIRECTORY, -out_dir OUTPUT_DIRECTORY
                        Enter the name of the directory where all other operations will be performed
  --genome GENOME, -g GENOME
                        Enter the SOFT-MASKED genome file of the organism

Optional arguments:
  --cpu CPU, -n CPU     Enter the number of CPUs to be used.
  --genome_dir_star GENOME_DIR_STAR, -gdir_star GENOME_DIR_STAR
                        Please enter the location of the genome index directory of STAR
  --genome_dir_olego GENOME_DIR_OLEGO, -gdir_olego GENOME_DIR_OLEGO
                        Please enter the location of the genome index directory of OLego
  --verbose VERBOSE, -verb VERBOSE
                        Enter a verbosity level
  --protein PROTEIN, -p PROTEIN
                        Enter the protein fasta
  --no_cleanup, -no_cleanup
                        Provide this option if you do not wish to remove any intermediate files. Please note that this will NOT remove any files and might take up a large amount of space
  --preserve_raw_input_data, -preserve
                        Set this argument if you want to preserve the raw fastq files. All other temporary files will be removed. These fastq files can be later used. 
  --checkpoint CHECKPOINT, -c CHECKPOINT
                        Enter a value if you wish to restart operations from a certain check point. Please note if you have new RNA-Seq samples, then FINDER will override this argument and computation will take place from read alignment. If there are missing data in any step then also FINDER will enforce restart of operations from a previous checkpoint. For example, if you wish to run assembly on samples for which alignments are not available then FINDER will readjust this value and set it to 1.
                            1. Align reads to reference genome (Will trigger removal of all alignments and start from beginning)
                            2. Assemble with PsiCLASS (Will remove all assemblies) 
                            3. Find genes with FINDER (entails changepoint detection) 
                            4. Predict genes using BRAKER2 (Will remove previous results of gene predictions with BRAKER2)
                            5. Annotate coding regions
                            6. Merge FINDER annotations with BRAKER2 predictions and protein sequences
```

FINDER can be launched using the following command:

```bash
finder -no_cleanup -mf Arabidopsis_thaliana_metadata.csv -n $CPU -gdir_star $PWD/star_index_without_transcriptome -out_dir $PWD/FINDER_test_ARATH -g $PWD/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa -p $PWD/uniprot_ARATH.fasta -gdir_olego olego_index -preserve 1> $PWD/FINDER_test_ARATH.output 2> $PWD/FINDER_test_ARATH.error 
```

This program will download and run the entire process of annotation. The duration of execution will depend on your internet speed and the number of cores you assigned to FINDER. Also, FINDER is designed in a way to handle a large number of RNA-Seq samples. So the speedup might not be noticeable with just a few samples.

### Enforcing running of FINDER from preset checkpoints

FINDER allows users to enforce execution from a specific checkpoints. Requesting a particular checkpoint does not mean that FINDER will skip all previous steps. It means that FINDER will remove all files generated by process after the checkpoint to ensure that the modules recalculate those. Below is a description of all the checkpoints that FINDER can accept:

1. Align reads to reference genome - Requesting FINDER to start from this checkpoint will trigger removal of all previous alignments. 
2. Assemble with PsiCLASS - Requesting FINDER to start from this checkpoint will trigger removal of assemblies that was previously generated. Aligned files will not be removed. If there are some RNA-Seq samples that are not aligned FINDER will align those first before attempting to assemble them
3. Find genes with FINDER - FINDER will regenerate all files post assembly by PsiCLASS
4. Predict genes using BRAKER2 - FINDER will rerun the BRAKER2 step
5. Annotate coding regions - FINDER will restart from annotating the coding sequences
6. Merge FINDER annotations with BRAKER2 predictions and protein sequences - FINDER will generate merged annotations from RNA-Seq samples, predictions and protein sequences

If you wish to start FINDER from downloading the SRA samples, please delete the output directory and start over.



## Output Files

All relevant output files generated by FINDER can be found in the `final_GTF_files` directory under the output directory. Below is the list of files and what data they contain

1. **braker.gtf** - gene models generated by BRAKER2
2. **combined_redundant_transcripts_removed.gtf** - GTF file from PsiCLASS output
3. **combined_split_transcripts_with_bad_SJ_redundancy_removed.gtf** - GTF file after splitting transcripts. This file is generated only from RNA-Seq expression evidence
4. **combined_with_CDS.gtf** - FINDER output *with* CDS predicted by GeneMark-S/T
5. **combined_with_CDS_high_conf.gtf**- FINDER gene models with high confidence
6. **combined_with_CDS_low_conf.gtf**- FINDER gene models with low confidence
7. **combined_with_CDS_BRAKER_appended_high_conf.gtf** - *High* confidence gene models from RNA-Seq evidence combined with BRAKER2 gene models
8. **combined_with_CDS_high_and_low_confidence_merged.gtf** - *High* and *Low* confidence gene models from RNA-Seq evidence combined with BRAKER2 gene models
9. **FINDER_BRAKER_PROT.gtf** - *High* confidence gene models from RNA-Seq evidence combined with BRAKER2 gene models and gene models from protein evidence
10. tissue/condition to transcript

## Intermediate files and folders

FINDER generates several intermediate files and folders. This section contains a detailed outline of the contents of each folder and what each file represents. 



## Utilities included with FINDER



## Terms of use



## FAQ



## Support

Create a [GitHub issue](https://github.com/sagnikbanerjee15/finder/issues).

