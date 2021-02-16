# Welcome to `finder`

`finder` is a gene annotator pipeline which automates the process of downloading short reads, aligning them and using the assembled  transcripts to generate gene annotations. Additionally it uses protein sequences and reports gene predictions by `BRAKER2`. It is a fast, scalable, platform independent software that generatess gene annotations in GTF format. `finder` accepts inputs through command line interface. It finds several novel genes/transcripts and also reports the tissue/conditions they were found to be in. If you use `finder` for your research please cite 

Sagnik Banerjee, Priyanka Bhandary, Margaret Woodhouse, Taner Z Sen,Roger P Wise, and Carson M Andorf.  FINDER: An automated software package to annotate eukaryotic genes from RNA-Seq data and associated protein sequences.bioRxiv, page 2021.02.04.429837, 2 2021.

## Installation

`finder` requires a number of softwares which needs to be installed. This might cause version conflicts with softwares that are already installed in your system. Hence, the developers have decided to enforce the use of `finder` within a conda environment. 

### Downloading and installing Anaconda (skip if you already have Anaconda installed)

Execute the following commands to download and install Anaconda

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
bash Anaconda3-2020.11-Linux-x86_64.sh # do not install VS. You may replace it with the latest version of Anaconda from https://www.anaconda.com/distribution/. Also, conda will default to the home directory, but make sure you choose a directory that has sufficient disk space to install all the software packages. 
```

## Downloading `finder` and creating a conda environment

```bash
git clone https://github.com/sagnikbanerjee15/finder.git
cd Finder
conda env create -f environment.yml # This will create an environment named finder_conda_env
conda activate finder_conda_env
cd dep
```

If you want to avoid using `git clone`and decide to download the zip package then please follow the following steps

```bash
wget https://github.com/sagnikbanerjee15/Finder/archive/master.zip
unzip master.zip
mv Finder-master Finder
cd Finder
conda env create -f environment.yml # This will create an environment named finder_conda_env
conda activate finder_conda_env
cd dep
```

`finder` runs `BRAKER2` which depends on `GeneMark-ES`. `finder` also needs `GeneMarkS/T` to predict coding regions of genes. Both `GeneMark-ES` and `GeneMarkS/T` are hosted at the University of Georgia website. The license prohibits the redistribution of their software, which is why it could not be included in this package. Hence, users have to manually download these 2 softwares and place them under the `dep` sub-directory. Please follow the instructions below to download the softwares and the key:

1. Open a browser of your choice
2. Go to [this](http://topaz.gatech.edu/GeneMark/license_download.cgi) website
3. Select the option **GeneMark-ES/ET/EP ver 4.62_lic** (2<sup>nd</sup> from top) and **LINUX 64**
4. Enter your name, institution, country and email-id and click on the button that says ***I agree to the terms of this license agreement***
5. Right click on the link that says *Please download program **here*** and select **Copy Link Address**
6. Then type in `wget ` and paste the path you just copied
7. This command will download the file **gmes_linux_64.tar.gz** in the `dep` sub-directory
8. Go to step 2
9. Select the option **GeneMarkS-T** (last option) and **LINUX 64**
10. Repeat steps 4-6
11. This command will download the file **gmst_linux_64.tar.gz** in the `dep` sub-directory
12. Now, right click on the link that says **64_bit** and select **Copy Link Address**
13. Then type in `wget` and paste the path you just copied
14. This command will download the file **gm_key_64.tar.gz** in the `dep` sub-directory. This license key will serve for both the programs. Please note that this key will expire after one year from the date of download.
15. Execute the following commands:

```bash
cd ..
./install.py
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc # Add this path permanently to the bashrc file
export PATH=$PATH:$(pwd)
```

## Executing FINDER with Sample data

Please follow the following the instructions to generate gene annotations using *Arabidopsis thaliana*. A `csv` file template has been provided with the release in `example/Arabidopsis_thaliana_metadata.csv`. Keep all the headers intact and replace the data with your samples of choice. Also note, that FINDER can work with both data downloaded from `NCBI` and also with data on local directories. Below is a detailed description of the each column of the metadata file.

| Column Name      | Column Description                                           | Mandatory |
| :--------------- | :----------------------------------------------------------- | :-------- |
| BioProject       | Name of the bioproject that the data belongs to. If you are using locally saved data then please enter a dummy project name. Please note that FINDER will **NOT** be able to process empty fields of Bioproject. | **YES**   |
| SRA Accession    | Enter the SRA Accession number of the sample that you expect FINDER to use for generating the gene annotations. Note that FINDER will use this ID to download the read samples from NCBI-SRA. In case you wish to use data which is not currently uploaded to NCBI, then you should enter the name of the local file. Do not enter any file extension in this field. For example, if your filename is `sample1.fastq`, please enter `sample1` in this field. `finder` assumes all files have the extension fastq. If there are files in your system that end with `f.q` please rename those to `*.fastq`. For paired-ended samples do not include the pair information in this field. For example, if you have 2 files `sample2_1.fastq` and `sample2_2.fastq` please enter `sample2` in this field. | **YES**   |
| Tissues          | Mention the tissue type or condition from which the sample has been collected. FINDER will report the tissues that are associated with a particular transcript. This can be used to find gene models that are expressed in a specific tissue and/or condition | **YES**   |
| Description      | A brief description of the data. This field is not mandatory and is not used by FINDER. It is upto the user to enter whatever metadata is deemed important. | NO        |
| Date             | Enter the date of producing the RNA-Seq sample. This field is not mandatory and is not used by FINDER. | NO        |
| Read Length (bp) | Enter the length of the reads. This field is not mandatory and is not used by FINDER. | NO        |
| Ended            | Enter either PE or SE for Paired ended reads or single neded reads. No other value should be entered. | **YES**   |
| RNA-Seq          | Enter 1 for all the rows. This field is included for future extensions. | **YES**   |
| process          | Enter 1 if you wish to process the sample. If a value of 0 is present, then FINDER will ignore the sample | **YES**   |
| Location         | Enter the location of the directory. For samples to be downloaded from NCBI, this field should be left empty. If the location of a directory is provided here then `finder` will assume that the sample is present in it. `finder` will generate an error if the sample is not found in this directory. It is *not* necessary to have all the samples in the same directory. | **YES**   |

To optimize disk space usage `finder` will process read samples from each bioproject at a time. Once the data is downloaded and reads are mapped, FINDER will remove all those data (if `-no-cleanup` is not specificied) to save disk space. But samples that were locally present will not be removed.

### Preparing the genome index

`finder` uses `STAR` and `OLego` for aligning and `PsiCLASS` for assembling. Users have the choice of generating the genome index and providing it to `finder` or allowing `finder` to generate the index. For some organisms with large genomes, `STAR` might require more memory. Hence users might be forced to transfer the STAR genome index generated in a computer that permits the usage of a large enough memory for index generation.  Download the *Arabidopsis thaliana* genome using the following command:

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
finder -h

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

`finder` can be launched using the following command:

```bash
finder -no_cleanup -mf Arabidopsis_thaliana_metadata.csv -n $CPU -gdir_star $PWD/star_index_without_transcriptome -out_dir $PWD/FINDER_test_ARATH -g $PWD/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa -p $PWD/uniprot_ARATH.fasta -gdir_olego olego_index -preserve 1> $PWD/FINDER_test_ARATH.output 2> $PWD/FINDER_test_ARATH.error 
```

This program will download and run the entire process of annotation. The duration of execution will depend on your internet speed and the number of cores you assigned to FINDER. Also, FINDER is designed in a way to handle a large number of RNA-Seq samples. So the speedup might not be noticeable with just a few samples.

Run the following command to remove all intermediate files. We recommend that while you run `finder`, you preserve all intermediate files and then run the following command to remove all the intermediate files.

```bash
finder -no_cleanup -mf Arabidopsis_thaliana_metadata.csv -n $CPU -gdir_star $PWD/star_index_without_transcriptome -out_dir $PWD/FINDER_test_ARATH -g $PWD/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa -p $PWD/uniprot_ARATH.fasta -gdir_olego olego_index -preserve -pc_clean 1> $PWD/FINDER_test_ARATH.output 2> $PWD/FINDER_test_ARATH.error
```

### Enforcing running of FINDER from preset checkpoints

`finder` allows users to enforce execution from a specific checkpoints. Requesting a particular checkpoint does not mean that `finder` will skip all previous steps. It means that `finder` will remove all files generated by process after the checkpoint to ensure that the modules recalculate those. Below is a description of all the checkpoints that `finder` can accept:

1. Align reads to reference genome - Requesting `finder` to start from this checkpoint will trigger removal of all previous alignments. 
2. Assemble with `PsiCLASS` - Requesting `finder` to start from this checkpoint will trigger removal of assemblies that was previously generated. Aligned files will not be removed. If there are some RNA-Seq samples that are not aligned FINDER will align those first before attempting to assemble them
3. Find genes with `finder` - `finder` will regenerate all files post assembly by PsiCLASS
4. Predict genes using `BRAKER2` - `finder` will rerun the BRAKER2 step
5. Annotate coding regions - `finder` will restart from annotating the coding sequences
6. Merge `finder` annotations with `BRAKER2` predictions and protein sequences - `finder` will generate merged annotations from RNA-Seq samples, predictions and protein sequences

If you wish to start `finder` from downloading the SRA samples, please delete the output directory and start over.

## Output Files

All relevant output files generated by `finder` can be found in the `final_GTF_files` directory under the output directory. Below is the list of files and what data they contain

1. **braker.gtf** - gene models generated by `BRAKER2`
2. **braker_utr.gtf** - gene models, with UTR models, generated by `BRAKER2`
3. **combined_redundant_transcripts_removed.gtf** - GTF file from PsiCLASS output
4. **combined_split_transcripts_with_bad_SJ_redundancy_removed.gtf** - GTF file after splitting transcripts. This file is generated only from RNA-Seq expression evidence
5. **combined_with_CDS.gtf** - `finder` output *with* CDS predicted by GeneMark-S/T
6. **combined_with_CDS_high_conf.gtf**- `finder` gene models with high confidence
7. **combined_with_CDS_low_conf.gtf**- `finder` gene models with low confidence
8. **combined_with_CDS_BRAKER_appended_high_conf.gtf** - *High* confidence gene models from RNA-Seq evidence combined with BRAKER2 gene models
9. **combined_with_CDS_high_and_low_confidence_merged.gtf** - *High* and *Low* confidence gene models from RNA-Seq evidence combined with BRAKER2 gene models
10. **FINDER_BRAKER_PROT.gtf** - *High* confidence gene models from RNA-Seq evidence combined with BRAKER2 gene models and gene models from protein evidence
11. **tissue/condition to transcript** - 

## Intermediate files and folders [To be updated]

`finder` generates several intermediate files and folders. This section contains a detailed outline of the contents of each folder and what each file represents. 

## Checking Progress

`finder` is configured to output information to a log file location in the output directory named `progress.log`. While reporting issues please make sure you attach the log file. 

## Restarting previous runs with more RNA-Seq samples

`finder` offeres users the opportunity to augment data into already completed annotation runs. Users need to update the metadata.csv file with the new RNA-Seq data and rerun `finder`. The program will determine an optimal starting point. `finder` will skip downloading of already processed RNA-Seq samples and will proceed with the new data. Users also have the option of removing some previously supplied RNA-Seq samples.

## Utilities included with FINDER

`finder` offers users with 2 utilites which could be used independently. 

1. `downloadAndDumpFastqFromSRA.py` - A python program that optimizes the download of data from SRA. Ids of RNA-Seq (or any sequencing for that matter) needs to be provided as a newline separated file. The program will download the RNA-Seq files, using the requested number of cores, convert those to fastq and remove the `.sra` files. `downloadAndDumpFastqFromSRA.py` will continuosly query the SRA database in the event of a failure.

   ```
   python downloadAndDumpFastqFromSRA.py -h
   usage: download_and_dump_fastq_from_SRA.py [-h] --sra SRA --output OUTPUT
                                              [--cpu CPU]
   
   Parallel download of fastq data from NCBI. Program will create the output
   directory if it is not present. If fastq file is present, then downloading is
   skipped. Program optimizes downloading of sra files and converting to fastq by
   utilizing multiple CPU cores.
   
   optional arguments:
     -h, --help            show this help message and exit
     --sra SRA, -s SRA     Please enter the name of the file which has all the
                           SRA ids listed one per line. Please note the
                           bioproject IDS cannot be processed
     --output OUTPUT, -o OUTPUT
                           Please enter the name of the output directory.
                           Download will be skipped if file is present
     --cpu CPU, -n CPU     Enter the number of CPUs to be used.
   ```

2. `verifyInputsToFINDER.py` - This program will verify whether all the resuested samples are in fact from a transcriptomic source of the organism whose genome is being annotated.

   ```
   python verifyInputsToFINDER.py -h
   usage: verify_inputs_to_finder.py [-h] --metadatafile METADATAFILE --srametadb
                                     SRAMETADB --taxon_id TAXON_ID
   
   Verifies whether all the data are transcriptomic and from the organism under
   consideration
   
   optional arguments:
     -h, --help            show this help message and exit
     --metadatafile METADATAFILE, -mf METADATAFILE
                           Please enter the name of the metadata file. Enter 0 in
                           the last column of those samples which you wish to
                           skip processing. The columns should represent the
                           following in order --> BioProject,Run,tissue_group,tis
                           sue,description,Date,read_length,ended (PE or
                           SE),directorypath,download,skip. If the sample is
                           skipped it will not be downloaded. Leave the directory
                           path blank if you are downloading the samples. In the
                           end of the run the program will output a csv file with
                           the directory path filled out. Please check the
                           provided csv file for more information on how to
                           configure the metadata file.
     --srametadb SRAMETADB, -m SRAMETADB
                           Enter the location of the SRAmetadb file.
     --taxon_id TAXON_ID, -t TAXON_ID
                           Enter the taxonomic id of the organism. Enter -1 if
                           you are working on a non-model organism or a sub-
                           species for which no taxonomic id exists.
   ```

   

## Terms of use

MIT License



Copyright (c) [2021] [Banerjee]



Permission is hereby granted, free of charge, to any person obtaining a copy

of this software and associated documentation files (the "Software"), to deal

in the Software without restriction, including without limitation the rights

to use, copy, modify, merge, publish, distribute, sublicense, and/or sell

copies of the Software, and to permit persons to whom the Software is

furnished to do so, subject to the following conditions:



The above copyright notice and this permission notice shall be included in all

copies or substantial portions of the Software.



THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR

IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,

FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE

AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER

LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,

OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE

SOFTWARE.

## Support

Please report all issues [here](https://github.com/sagnikbanerjee15/finder/issues).

