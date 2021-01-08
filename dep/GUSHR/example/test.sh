#!/bin/bash

if [ ! -f RNAseq.bam ]; then
    wget http://bioinf.uni-greifswald.de/bioinf/braker/RNAseq.bam
    fi

../gushr.py -b RNAseq.bam -t augustus.gtf -g genome.fa -o gushr -c 10
