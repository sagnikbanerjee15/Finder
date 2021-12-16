

# Change Log

All notable changes to this project will be documented here. This file is divided into different groups. The unreleased section contains updates that are currently being processed in the lab. The future enhancement section contains updates that have been requested by users and/or have been discussed in the lab but is not being currently developed.

## [1.1.0] - Jan 10<sup>th</sup>, 2022

Here we write upgrading notes for `finder`. It's an effort to make them as straightforward as possible.

### Added

- Option to provide alignments by `exonerate` to `finder`
- `finder` will now be able to process gzipped `fastq` files. But this option will work only if the files are available locally

### Changed

* Execution within`docker` and `singularity` options incorporated
* `GeneMark` software path and license has to be provided as arguments to the main program
* `finder` will generate `STAR` and `Olego` indices
* Coding Sequence annotation was improved by incorporating valid ORFs
* `environment.yml` file was removed since execution in conda environment is no longer being offered
* Removed `setup.py` file. 
* Total time for downloading, aligning and assembling all samples provided in `progress.log` file
* Coding sequence prediction software has been changed from `GeneMarkS/T` to `CodAn`

### Fixed

* Issues with reading the `metadata.csv` file, especially when some fields are blank, have been resolved

## Future enhancements (FINDER2)

* Add option to predict transcription factor binding sites from motif data by incorporating softwares from. MEME suite
* Add option to process different kinds of NGS data like CAGE-Seq, RAMPAGE-Seq, Ribo-Seq etc
* Modify the CPD package. Also recode the parts in C programming language to make things run fast
* Try gene predictors other than `GeneMark` to circumvent license issue
* Also change CDS prediction technique
* Simpler alignment strategy - do away with 4 rounds of mapping. This will increase speed considerably. 
* Replace OLego with some de-novo strategy of assembling. Olego is currenlty taking way too much time.



## [1.0.0] - July 2<sup>nd</sup>, 2021

* First release of `finder`
* Please download it from https://github.com/sagnikbanerjee15/Finder/releases/tag/finder_v1.0.0

