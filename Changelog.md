

# Change Log

All notable changes to this project will be documented here. This file is divided into different groups. The unreleased section contains updates that are currently being processed in the lab. The future enhancement section contains updates that have been requested by users and/or have been discussed in the lab but is not being currently developed.

## [Unreleased] - 2021-03-02

Here we write upgrading notes for `finder`. It's an effort to make them as straightforward as possible.

### Added

- Option to provide alignments by `exonerate` to `finder`
- `finder` will now be able to process gzipped `fastq` files. But this option will work only if the files are available locally
- `finder` performs enhanced memory checks to ensure 
- Utility program added to convert `gtf` files to `gff3` including options to generate UTR annotation. `README.md` file has been updated accordingly
- Included an option to use local RNA-Seq files in the testing data

### Changed

* Coding Sequence annotation was improved by incorporating valid ORFs
* `finder` now verifies the maximum length of command and issues a warning
* `psiclass` developers were requested to modify the program to incorporate options to adjust the length of the end-exons based on RNA-Seq coverage data. The latest version of `psiclass` is available from conda 
* Updated `environment.yml` file. `psiclass` will now be installed directly from `conda`
* Updated `setup.py` file. Removed the part where `psiclass` was being installed
* Option added to skip `braker` run completely
* Option added to incorporate reads from PacBio or any other long read technology
* Functionality added to merge genes that are close to one another (separated by only a few nucleotides). Some gene models of `psiclass` have been split despite a continuous coverage
* Add options to perform repeat masking of unmasked genome
* Total time for downloading, aligning and assembling all samples provided in `progress.log` file

### Fixed

* Issues with reading the `metadata.csv` file, especially when some fields are blank, have been resolved

### Deprecated



### Removed



### Security



## Future enhancements (FINDER2)

* Add option to predict transcription factor binding sites from motif data by incorporating softwares from. MEME suite
* Add option to process different kinds of NGS data like CAGE-Seq, RAMPAGE-Seq, Ribo-Seq etc
* Modify the CPD package. Also recode the parts in C programming language to make things run fast
* Try gene predictors other than `GeneMark` to circumvent license issue
* Also change CDS prediction technique
* Simpler alignment strategy - do away with 4 rounds of mapping. This will increase speed considerably. 
* Replace OLego with some de-novo strategy of assembling. Olego is currenlty taking way too much time.



## [1.0.0] - 2021-02-06

* First release of `finder`

