#! /usr/bin/env python3

# Quick and dirty - Need to change later

import copy
import math
import os
import sys

gene_annotation_gtf = sys.argv[1]

# Read from gtf file

fhr = open(gene_annotation_gtf,"r")
gene_info = {}
for line in fhr:
    line = line.strip().split("\t")
    if line[2] == "transcript":
        chromosome = line[0]
        if chromosome not in gene_info:
            gene_info[chromosome] = []
        previous_transcript_id = line[-1].split("transcript_id")[-1].strip().split()[0].strip("\"")
        start = int(line[3])
        end = int(line[4])
        strand = line[6]
        gene_info[chromosome].append([previous_transcript_id,start,end,strand])
fhr.close()

# Merge transcripts into genes

transcripts_grouped_into_genes = {}

for chromosome in gene_info:
    gene_number = 1
    for row in gene_info[chromosome]:
        previous_transcript_id, start, end, strand = row
        transcript_incorporated = 0
        if chromosome not in transcripts_grouped_into_genes:
            transcripts_grouped_into_genes[chromosome] = {}
        for gene in transcripts_grouped_into_genes[chromosome]:
            for transcript in transcripts_grouped_into_genes[chromosome][gene]:
                transcript_start, transcript_end, transcript_strand = transcripts_grouped_into_genes[chromosome][gene][transcript]
                if strand !="." and strand != transcript_strand: continue 
                # Check if this transcript overlaps at least 80%
                if start <= transcript_start <= end <= transcript_end:
                    if ( end - transcript_start )/(end - start if (end-start)<(transcript_end-transcript_start) else (transcript_end-transcript_start) ) > 0.8:
                        transcripts_grouped_into_genes[chromosome][gene].append(previous_transcript_id)
                        transcript_incorporated = 1
                elif transcript_start <= start <= transcript_end <= end:
                    if ( transcript_end - start )/(end - start if (end-start)<(transcript_end-transcript_start) else (transcript_end-transcript_start) ) > 0.8:
                        transcripts_grouped_into_genes[chromosome][gene].append(previous_transcript_id)
                        transcript_incorporated = 1
                elif start <= transcript_start <= transcript_end <= end:
                    transcripts_grouped_into_genes[chromosome][gene].append(previous_transcript_id)
                    transcript_incorporated = 1
                elif transcript_start <= start <= end <= transcript_end:
                    transcripts_grouped_into_genes[chromosome][gene].append(previous_transcript_id)
                    transcript_incorporated = 1
        if transcript_incorporated == 0:
            transcripts_grouped_into_genes[chromosome][f"gene_{gene_number}"] = [previous_transcript_id]
            gene_number+=1
            
for chromosome in transcripts_grouped_into_genes:
    for gene in transcripts_grouped_into_genes[chromosome]:
        print(chromosome, gene, transcripts_grouped_into_genes[chromosome][gene])