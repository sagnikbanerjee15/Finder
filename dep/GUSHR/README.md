# GUSHR

Assembly-free construction of UTRs from short read RNA-Seq data on the basis of coding sequence annotation.

This tool has been adapted to the format needs of AUGUSTUS/BRAKER and employs GeMoMa for generating UTRs from RNA-Seq coverage data.

## Contacts for Github Repository of GUSHR:

Katharina J. Hoff, University of Greifswald, Germany, katharina.hoff@uni-greifswald.de, +49 3834 420 4624

## Authors of GUSHR:

Katharina J. Hoff<sup name="aff1">[a, ](#aff1)</sup><sup name="aff2">[b](#aff2)</sup>, Mario Stanke <sup name="aff1">[a, ](#aff1)</sup><sup name="aff2">[b](#aff2)</sup> and Jens Keilwagen<sup name="aff3">[c](#aff3)</sup>

<b id="aff1">[a]</b> University of Greifswald, Institute for Mathematics and Computer Science, Walther-Rathenau-Str. 47, 17489 Greifswald, Germany

<b id="aff2">[b]</b> University of Greifswald, Center for Functional Genomics of Microbes, Felix-Hausdorff-Str. 8, 17489 Greifswald, Germany

<b id="aff3">[c]</b> Julius Kühn-Institut, Erwin-Baur-Str. 27, 06484 Quedlinburg, Germany


## Software Requirements

GUSHR is a script written in Python3 that runs a number of external tools, some of which are Bash core utils. Thus, GUSHR should be executed on a Linux system (we developed and tested code on Ubuntu 20.04.1 LTS).

Required Bash core utils:

   * grep

   * sort

Other software dependencies:

   * java 1.8

   * samtools 1.8-20-g4ff8062

   * gtf2gff.pl from AUGUSTUS 3.3.3 (https://github.com/Gaius-Augustus/Augustus)

GUSHR runs a GeMoMa.jar file that is supplied with GUSHR. This GeMoMa.jar file was originally created by Jens Keilwagen et al. with Java 1.8 (http://www.jstacs.de/index.php/GeMoMa#Download, https://github.com/Jstacs/Jstacs) and is here re-distributed. You have the option to specify a different (e.g. more recent) GeMoMa.jar file when calling GUSHR but compatibility is then not guaranteed.

## Test GUSHR

First, you need to retrieve an RNA-Seq example file that is not included on github due to its size:

    wget http://bioinf.uni-greifswald.de/bioinf/braker/RNAseq.bam

(or: `cd example; ./download_rnaseq.sh`).

After that, and if all software requirements are satisfied, you can test GUSHR:

    ./gushr.py -t example/augustus.gtf -b example/RNAseq.bam -g example/genome.fa -o utrs

## What to Cite

If you used GUSHR in your work, please cite the following sources:

   * Stanke, M., Diekhans, M., Baertsch, R. and Haussler, D. (2008). Using native and syntenically mapped cDNA alignments to improve de novo gene finding. Bioinformatics, doi: 10.1093/bioinformatics/btn013

   * Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R.; 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16):2078-9.

   * Hoff, K. J., & Stanke, M. (2019). Predicting genes in single genomes with augustus. Current protocols in bioinformatics, 65(1), e57.

   * Keilwagen, J., Hartung, F., Grau, J. (2019) GeMoMa: Homology-based gene prediction utilizing intron position conservation and RNA-seq data. Methods Mol Biol. 1962:161-177, doi: 10.1007/978-1-4939-9173-0_9.

   * Keilwagen, J., Wenk, M., Erickson, J.L., Schattat, M.H., Grau, J., Hartung F. (2016) Using intron position conservation for homology-based gene prediction. Nucleic Acids Research, 44(9):e89.

   * Keilwagen, J., Hartung, F., Paulini, M., Twardziok, S.O., Grau, J. (2018) Combining RNA-seq data and homology-based gene prediction for plants, animals and fungi. BMC Bioinformatics, 19(1):189.