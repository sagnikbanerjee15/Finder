#! /usr/bin/env python3

# Quick and dirty - Need to change later

import copy
import math
import os
import sys

import pandas as pd

gene_annotation_gtf = sys.argv[1]

# Read from gtf file

fhr = open(gene_annotation_gtf,"r")
gene_info = {}
for line in fhr:
    line = line.strip().split("\t")
    if line[2] == "transcript":
        chromosome = line[0]
        if chromosome not in gene_info:
            gene_info[chromosome] = {}
        previous_transcript_id = line[-1].split("transcript_id")[-1].strip().split()[0].strip("\"")
        start = int(line[3])
        end = int(line[4])
        gene_info[chromosome][previous_transcript_id] = [start,end]
fhr.close()