#! /usr/bin/env python

from scripts.fileReadWriteOperations import *
import copy
import math
import os
import sys

import pandas as pd


def mergeTwoTranscripts(whole_annotations,transcript_id_i,transcript_id_j,chromosome):
    """
    """
    #print("Merging",transcript_id_i,transcript_id_j)
    chromosome=transcript_id_i.split(".")[0]
    transcript_id_i_info=whole_annotations[transcript_id_i]
    transcript_id_j_info=whole_annotations[transcript_id_j]
    new_transcript_id=".".join(transcript_id_i.split(".")[:-1])+"_"+transcript_id_i.split(".")[-1]+"_merged_"+"_".join(transcript_id_j.split(".")[:-1])+"."+transcript_id_j.split(".")[-1]
    #print(transcript_id_i,transcript_id_j,new_transcript_id)
    sys.stdout.flush()
    whole_annotations[new_transcript_id]={"exons":copy.deepcopy(whole_annotations[transcript_id_i]["exons"]),
                                                      "introns":[],
                                                      "cov":whole_annotations[transcript_id_i]["cov"],
                                                      "TPM":whole_annotations[transcript_id_i]["TPM"],
                                                      "FPKM":whole_annotations[transcript_id_i]["FPKM"],
                                                      "direction":whole_annotations[transcript_id_i]["direction"],
                                                      "chromosome":chromosome
                                                      }
    
    whole_annotations[new_transcript_id]["exons"][-1]=[whole_annotations[transcript_id_i]["exons"][-1][0],
                                                                   whole_annotations[transcript_id_j]["exons"][0][1]]
    
    if len(whole_annotations[transcript_id_j]["exons"])>1:
        whole_annotations[new_transcript_id]["exons"].extend(whole_annotations[transcript_id_j]["exons"][1:])
    
    i=1
    while i<len(whole_annotations[new_transcript_id]["exons"]):
        whole_annotations[new_transcript_id]["introns"].append([whole_annotations[new_transcript_id]["exons"][i-1][1]+1,whole_annotations[new_transcript_id]["exons"][i][0]-1])
        i+=1
    return whole_annotations
 
 
def mergeCloselySpacedTranscripts(options):
    """
    """
    input_gtf_filename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_cov_opp_split_redundancy_removed.gtf"
    output_gtf_filename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_merged_transcripts.gtf"
    if os.path.exists(output_gtf_filename)==True:return
    whole_annotations,useless1,useless2=readAllTranscriptsFromGTFFileInParallel([input_gtf_filename,"dummy","dummy"])
    all_transcript_info=[]
    
    for transcript_id in whole_annotations:
        chromosome = whole_annotations[transcript_id]["chromosome"]
        transcript_start=whole_annotations[transcript_id]["transcript_start"]
        transcript_end=whole_annotations[transcript_id]["transcript_end"]
        cov=whole_annotations[transcript_id]["cov"]
        fpkm=whole_annotations[transcript_id]["FPKM"]
        tpm=whole_annotations[transcript_id]["TPM"]
        direction=whole_annotations[transcript_id]["direction"]
        all_transcript_info.append([chromosome,transcript_id,transcript_start,transcript_end,cov,fpkm,tpm,direction])
    all_transcript_info_pd=pd.DataFrame(all_transcript_info,columns=["chromosome","transcript_id","transcript_start","transcript_end","cov","fpkm","tpm","direction"])
    all_transcript_info_pd=all_transcript_info_pd.sort_values(by=["chromosome","transcript_start"])
    remove_these_transcripts=[]
    for row_num,row in all_transcript_info_pd.iterrows():
        chromosome,transcript_id,transcript_start,transcript_end,cov,fpkm,tpm,direction=row
        if direction==".":continue
        potential_merger_transcript=all_transcript_info_pd[(all_transcript_info_pd["chromosome"]==chromosome) &
                               (all_transcript_info_pd["transcript_id"]!=transcript_id) & 
                               (all_transcript_info_pd["transcript_start"]>=transcript_end) &
                               (all_transcript_info_pd["direction"]==direction) &
                               (all_transcript_info_pd["transcript_start"]-transcript_end<=5 )
                               ]
        if potential_merger_transcript.shape[0]>0:
            for row_num_i,row_i in potential_merger_transcript.iterrows():
                chromosome_i,transcript_id_i,transcript_start_i,transcript_end_i,cov_i,fpkm_i,tpm_i,direction_i=row_i
                if math.fabs(tpm-tpm_i)<2 and max(tpm,tpm_i)<5 and "cov" not in transcript_id and "cov" not in transcript_id_i:
                    #print(transcript_id,transcript_id_i,tpm,tpm_i)
                    remove_these_transcripts.append(transcript_id)
                    remove_these_transcripts.append(transcript_id_i)
                    whole_annotations=mergeTwoTranscripts(whole_annotations,transcript_id,transcript_id_i,chromosome_i)
                sys.stdout.flush()
    
    for transcript_id in list(set(remove_these_transcripts)):
        chromosome=transcript_id.split(".")[0]
        del whole_annotations[transcript_id]
    
    writeTranscriptsToFile([whole_annotations,output_gtf_filename,0])
    
    