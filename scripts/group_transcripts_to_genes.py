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
        previous_transcript_id = line[-1].split("transcript_id")[-1].strip().split()[0].strip()[:-1].strip("\"")
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
                        transcripts_grouped_into_genes[chromosome][gene][previous_transcript_id] = [start, end, strand]
                        transcript_incorporated = 1
                        break
                elif transcript_start <= start <= transcript_end <= end:
                    if ( transcript_end - start )/(end - start if (end-start)<(transcript_end-transcript_start) else (transcript_end-transcript_start) ) > 0.8:
                        transcripts_grouped_into_genes[chromosome][gene][previous_transcript_id] = [start, end, strand]
                        transcript_incorporated = 1
                        break
                elif start <= transcript_start <= transcript_end <= end:
                    transcripts_grouped_into_genes[chromosome][gene][previous_transcript_id] = [start, end, strand]
                    transcript_incorporated = 1
                    break
                elif transcript_start <= start <= end <= transcript_end:
                    transcripts_grouped_into_genes[chromosome][gene][previous_transcript_id] = [start, end, strand]
                    transcript_incorporated = 1
                    break
        if transcript_incorporated == 0:
            transcripts_grouped_into_genes[chromosome][gene_number] = {}
            transcripts_grouped_into_genes[chromosome][gene_number][previous_transcript_id] = [start, end, strand]
            gene_number+=1

final_mapping_prev_transcript_to_new_gene_and_new_transcript_ids = {}
            
for chromosome in transcripts_grouped_into_genes:
    if chromosome not in final_mapping_prev_transcript_to_new_gene_and_new_transcript_ids:
        final_mapping_prev_transcript_to_new_gene_and_new_transcript_ids[chromosome] = {}
    for gene in transcripts_grouped_into_genes[chromosome]:
        for transcript_number,transcript in enumerate(transcripts_grouped_into_genes[chromosome][gene]):
            transcripts_grouped_into_genes[chromosome][gene][transcript] = transcript_number + 1
            final_mapping_prev_transcript_to_new_gene_and_new_transcript_ids[chromosome][transcript] = [gene, transcripts_grouped_into_genes[chromosome][gene][transcript]]
        
        #print(chromosome, gene, transcripts_grouped_into_genes[chromosome][gene])
        
for chromosome in final_mapping_prev_transcript_to_new_gene_and_new_transcript_ids:
    for transcript in final_mapping_prev_transcript_to_new_gene_and_new_transcript_ids[chromosome]:
        pass
        #print(chromosome, transcript, final_mapping_prev_transcript_to_new_gene_and_new_transcript_ids[chromosome][transcript])

fhw = open(gene_annotation_gtf[:-4] + "_final" + ".gtf", "w")
fhr = open(gene_annotation_gtf,"r")
for line in fhr:
    line = line.strip().split("\t")
    old_transcript_id = line[-1].split("transcript_id")[-1].strip().split()[0].strip()[:-1].strip("\"")
    old_gene_id = line[-1].split("gene_id")[-1].strip().split()[0].strip()[:-1].strip("\"")
    chromosome = line[0]
    new_gene_id, new_transcript_id = final_mapping_prev_transcript_to_new_gene_and_new_transcript_ids[chromosome][old_transcript_id]
    line[-1] = f"gene_id \"{chromosome}.{new_gene_id}\"; transcript_id \"{chromosome}.{new_gene_id}.{new_transcript_id}\"; "
    fhw.write("\t".join(line) + "\n")

fhr.close()
fhw.close()