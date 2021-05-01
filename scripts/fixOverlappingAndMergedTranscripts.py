

from copy import deepcopy
from scripts.fileReadWriteOperations import *
from scripts.runCommand import *
import multiprocessing
import pickle
import pprint
import sys

import numpy as np
import pandas as pd


def divide_chunks(l, n): 
      
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 

def removeSpuriousExonsAndTranscripts(whole_annotation):
    remove_these_transcripts = []
    for transcript_id in whole_annotation:
        # Remove transcripts that are less than 50 nucleotides long
        if whole_annotation[transcript_id]["transcript_end"] - whole_annotation[transcript_id]["transcript_start"] + 1 <50 :
            remove_these_transcripts.append(transcript_id)
        
    return whole_annotation
        
def createNewTranscripts(combined_list_of_exons_and_breakpoints,transcript_id):
    """
    """
    
    # Merge consequent cut-points
    combined_list_of_exons_and_breakpoints_new = []
    while True:
        i=0
        while i<len(combined_list_of_exons_and_breakpoints):
            if i+1<len(combined_list_of_exons_and_breakpoints) and combined_list_of_exons_and_breakpoints[i][-1]=='c' and combined_list_of_exons_and_breakpoints[i+1][-1]=='c':
                combined_list_of_exons_and_breakpoints_new.append([combined_list_of_exons_and_breakpoints[i][0],
                                                                   combined_list_of_exons_and_breakpoints[i+1][1],'c'])
                i+=1
            else:
                combined_list_of_exons_and_breakpoints_new.append(combined_list_of_exons_and_breakpoints[i])
            i+=1 
        """print(transcript_id,len(combined_list_of_exons_and_breakpoints),len(combined_list_of_exons_and_breakpoints_new))
        sys.stdout.flush()"""
        combined_list_of_exons_and_breakpoints = combined_list_of_exons_and_breakpoints_new
        i=0
        flag=0
        while i<len(combined_list_of_exons_and_breakpoints)-1:
            if  combined_list_of_exons_and_breakpoints[i][-1]=='c' and combined_list_of_exons_and_breakpoints[i+1][-1]=='c':
                flag=1
                break
            i+=1
        if flag==0:
            break
        combined_list_of_exons_and_breakpoints_new = []
        
    # Remove cut points that span both introns and exons
    combined_list_of_exons_and_breakpoints_new = []
    combined_list_of_exons_and_breakpoints_new.append(combined_list_of_exons_and_breakpoints[0])
    i=0
    while i<len(combined_list_of_exons_and_breakpoints):
        if combined_list_of_exons_and_breakpoints[i][-1]=='c':
            exon1_start,exon1_end = combined_list_of_exons_and_breakpoints[i-1][:2]
            cutpoint_start,cutpoint_end = combined_list_of_exons_and_breakpoints[i][:2]
            if exon1_start<=cutpoint_start and cutpoint_end<=exon1_end:
                combined_list_of_exons_and_breakpoints_new.append(combined_list_of_exons_and_breakpoints[i])
            else:
                pass
        else:
            combined_list_of_exons_and_breakpoints_new.append(combined_list_of_exons_and_breakpoints[i])
        i+=1
    combined_list_of_exons_and_breakpoints = combined_list_of_exons_and_breakpoints_new
    
    # Remove those cases where cut points fall at extreme edges (within 5 nucl) of exon boundaries
    combined_list_of_exons_and_breakpoints_new = []
    i=0
    while i<len(combined_list_of_exons_and_breakpoints):
        if combined_list_of_exons_and_breakpoints[i][-1]=='c':
            if i-1>=0 and i+1<len(combined_list_of_exons_and_breakpoints)-1:
                prev_exon = combined_list_of_exons_and_breakpoints[i-1][:2]
                next_exon = combined_list_of_exons_and_breakpoints[i+1][:2]
                cut_point_start, cut_point_end = combined_list_of_exons_and_breakpoints[i][:2]
                if abs(cut_point_start-prev_exon[0])<=5 or abs(cut_point_end-next_exon[1])<=5 or abs(cut_point_start-next_exon[0])<=5 or abs(cut_point_end-prev_exon[1])<=5:
                    pass
                else:
                    combined_list_of_exons_and_breakpoints_new.append(combined_list_of_exons_and_breakpoints[i])
            elif i == 0:
                next_exon = combined_list_of_exons_and_breakpoints[i+1][:2]
                cut_point_start, cut_point_end = combined_list_of_exons_and_breakpoints[i][:2]
                if abs(cut_point_end-next_exon[1])<=5 or abs(cut_point_start-next_exon[0])<=5:
                    pass
                else:
                    combined_list_of_exons_and_breakpoints_new.append(combined_list_of_exons_and_breakpoints[i])
            elif i == len(combined_list_of_exons_and_breakpoints)-1:
                prev_exon = combined_list_of_exons_and_breakpoints[i-1][:2]
                cut_point_start, cut_point_end = combined_list_of_exons_and_breakpoints[i][:2]
                if abs(cut_point_start-prev_exon[0])<=5 or abs(cut_point_end-prev_exon[1])<=5:
                    pass
                else:
                    combined_list_of_exons_and_breakpoints_new.append(combined_list_of_exons_and_breakpoints[i])
        else:
            combined_list_of_exons_and_breakpoints_new.append(combined_list_of_exons_and_breakpoints[i])
        i+=1
    #combined_list_of_exons_and_breakpoints_new.append( combined_list_of_exons_and_breakpoints[-1])
    """print("Prev_length",len(combined_list_of_exons_and_breakpoints))
    print("Next_length",len(combined_list_of_exons_and_breakpoints_new))"""
    combined_list_of_exons_and_breakpoints = combined_list_of_exons_and_breakpoints_new
    
    if len(combined_list_of_exons_and_breakpoints)>=2:
        combined_list_of_exons_and_breakpoints_new = []
        if combined_list_of_exons_and_breakpoints[-1][2]=='c' and combined_list_of_exons_and_breakpoints[-2][1] == combined_list_of_exons_and_breakpoints[-1][1] :
            combined_list_of_exons_and_breakpoints_new = combined_list_of_exons_and_breakpoints[:-1]
            combined_list_of_exons_and_breakpoints = combined_list_of_exons_and_breakpoints_new
    
    subtract_this_value = combined_list_of_exons_and_breakpoints[0][0]
    combined_list_of_exons_and_breakpoints_new = []
    for row in combined_list_of_exons_and_breakpoints:
        combined_list_of_exons_and_breakpoints_new.append([row[0]-subtract_this_value,
                                                           row[1]-subtract_this_value,
                                                           row[2]])
        
    combined_list_of_exons_and_breakpoints = combined_list_of_exons_and_breakpoints_new
    list_of_exon_definitons_of_transcripts=[]
    current_transcript_exon_definition=[]
    exon_num=0
    while exon_num<len(combined_list_of_exons_and_breakpoints):
        flag=0
        exon=combined_list_of_exons_and_breakpoints[exon_num]
        if exon[2]=='e' and exon_num+1<len(combined_list_of_exons_and_breakpoints) and combined_list_of_exons_and_breakpoints[exon_num+1][2]=='e':
            current_transcript_exon_definition.append([exon[0],exon[1]])
        elif exon[2]=='e' and exon_num+1<len(combined_list_of_exons_and_breakpoints) and combined_list_of_exons_and_breakpoints[exon_num+1][2]=='c':
            # Check if the cutpoint is overlapping with intron
            cutpoint = combined_list_of_exons_and_breakpoints[exon_num+1]
            exon1 = exon
            try:
                exon2 = combined_list_of_exons_and_breakpoints[exon_num+2]
            except IndexError:
                #print(transcript_id)
                #pprint.pprint(combined_list_of_exons_and_breakpoints)
                #sys.exit()
                flag=1
            if flag==0 and exon1[1]<=cutpoint[0] and cutpoint[1]<=exon2[0]:
                current_transcript_exon_definition.append([exon[0],exon[1]])
                current_transcript_exon_definition_cp=deepcopy(current_transcript_exon_definition)
                list_of_exon_definitons_of_transcripts.append(current_transcript_exon_definition_cp)
                current_transcript_exon_definition=[[exon2[0],exon2[1]]]
            else:
                current_transcript_exon_definition.append([exon[0],combined_list_of_exons_and_breakpoints[exon_num+1][0]-1])
                current_transcript_exon_definition_cp=deepcopy(current_transcript_exon_definition)
                list_of_exon_definitons_of_transcripts.append(current_transcript_exon_definition_cp)
                current_transcript_exon_definition=[[combined_list_of_exons_and_breakpoints[exon_num+1][1]+1,exon[1]]]
            exon_num+=1
        elif exon[2]=='e' and exon_num==len(combined_list_of_exons_and_breakpoints)-1:
            current_transcript_exon_definition.append([exon[0],exon[1]])
        exon_num+=1 
    #current_transcript_exon_definition.append(current_transcript_exon_definition)
    list_of_exon_definitons_of_transcripts.append(current_transcript_exon_definition)
    
    #print(transcript_id)
    """if transcript_id=="1.7600.0" or transcript_id=="3.41089.0" or transcript_id=="3.41089.1" or transcript_id=="3.46058.0" or transcript_id=="4.47251.1" or transcript_id=="4.47251.2"  or transcript_id=="5.71927.0"  or transcript_id=="5.71927.1":
        print(transcript_id)
        pprint.pprint(combined_list_of_exons_and_breakpoints)
        for each_transcript_exon_definition in list_of_exon_definitons_of_transcripts:
            pprint.pprint(each_transcript_exon_definition)
            print("-"*150)
        print("="*150)"""
     
    list_of_exon_definitons_of_transcripts_new=[]
    for each_transcript_exon_definition in list_of_exon_definitons_of_transcripts:
        each_transcript_exon_definition_new = []
        for row in each_transcript_exon_definition:
            each_transcript_exon_definition_new.append([row[0]+subtract_this_value,
                                                        row[1]+subtract_this_value])
        list_of_exon_definitons_of_transcripts_new.append(each_transcript_exon_definition_new)
    list_of_exon_definitons_of_transcripts = list_of_exon_definitons_of_transcripts_new
    
    return list_of_exon_definitons_of_transcripts

def removeOverlappingExonsFromEachTranscript(transcript_info):
    """
    Scans each transcript and removes those exons that are subsets of other exons
    """
    for transcript_id in transcript_info:
        remove_these_exons = []
        i=0
        while i<len(transcript_info[transcript_id]["exons"]):
            j=i+1
            exon_i = transcript_info[transcript_id]["exons"][i]
            while j<len(transcript_info[transcript_id]["exons"]):
                exon_j = transcript_info[transcript_id]["exons"][j]
                if exon_i[0]==exon_j[0] and exon_i[1]==exon_j[1]:
                    remove_these_exons.append(i)
                # Check if ith exon is contained in jth exon
                elif exon_j[0]<=exon_i[0] and exon_i[1]<=exon_j[1]: 
                    remove_these_exons.append(i)
                # Check if jth exon is contained in ith exon
                elif exon_i[0]<=exon_j[0] and exon_j[1]<=exon_i[1]: 
                    remove_these_exons.append(j)
                j+=1
            i+=1
        remove_these_exons = list(set(remove_these_exons))[::-1]
        """if transcript_id in ["1.15397_0_covsplit.0 ","1.17646_0_covsplit.0","1.19072_0_covsplit.0","4.50055_0_covsplit.0","5.68350_0_covsplit.0","5.68602_0_covsplit.0","5.69132_0_covsplit.0","5.71535_0_covsplit.0"]:
            print(transcript_id,len(transcript_info[transcript_id]["exons"]),len(transcript_info[transcript_id]["exons"])-len(remove_these_exons))
            print(transcript_info[transcript_id]["exons"])
            print("="*150)"""
        for r in remove_these_exons:
            transcript_info[transcript_id]["exons"].pop(r)
        
    return transcript_info
    
def fixOverlappingAndMergedTranscripts(options,logger_proxy,logging_mutex):
    
    #########################################################################################################
    # Read in all the exons that overlap with some intron
    #########################################################################################################
    
    # Selecting the first Run from the first condition
    gtffilename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_transcripts_connecting_two_transcripts.gtf"
    output_gtf_filename = options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_cov_opp_split.gtf"
    if options.skip_cpd==True:
        os.system(f"mv {gtffilename} {output_gtf_filename}")
        return
    if os.path.exists(output_gtf_filename)==True:return
    exons_overlapping_with_introns_bedfilename=gtffilename[:-4]+"_exons_overlapping_with_introns.bed"
    overlaps={}
    
    cmd="bedtools getfasta -s -fi "+options.genome+" -bed <(cat "+gtffilename+"|awk '$3==\"exon\"') > "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/sequence_per_exon.fasta"
    runCommand(["dummy",cmd])
    
    donot_split_these_exons={}
    sequence_per_exon=readFastaFile(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/sequence_per_exon.fasta")
    for id in sequence_per_exon:
        if len(sequence_per_exon[id])%3==0:
            frame_0=translate(sequence_per_exon[id])
        elif len(sequence_per_exon[id])%3==1:
            frame_0=translate(sequence_per_exon[id][:-1])
        elif len(sequence_per_exon[id])%3==2:
            frame_0=translate(sequence_per_exon[id][:-2])
        
        if len(sequence_per_exon[id][1:])%3==0:    
            frame_1=translate(sequence_per_exon[id][1:])
        elif len(sequence_per_exon[id][1:])%3==1:    
            frame_1=translate(sequence_per_exon[id][1:-1])
        elif len(sequence_per_exon[id][1:])%3==2:    
            frame_1=translate(sequence_per_exon[id][1:-2])
        
        if len(sequence_per_exon[id][2:])%3==0:
            frame_2=translate(sequence_per_exon[id][2:])
        elif len(sequence_per_exon[id][2:])%3==1:
            frame_2=translate(sequence_per_exon[id][2:-1])
        elif len(sequence_per_exon[id][2:])%3==2:
            frame_2=translate(sequence_per_exon[id][2:-2])
        
        skip=0
        if abs(len(sequence_per_exon[id])/3-len(frame_0))<=1:
            skip=1
        if abs(len(sequence_per_exon[id])/3-len(frame_1))<=1:
            skip=1
        if abs(len(sequence_per_exon[id])/3-len(frame_2))<=1:
            skip=1
        if skip==1:
            chromosome=id.strip()[:-3].split(":")[0]
            start=str(int(id.strip()[:-3].split(":")[-1].split("-")[0])+1)
            end=id.strip()[:-3].split(":")[-1].split("-")[-1]
            if chromosome not in donot_split_these_exons:
                donot_split_these_exons=[]
            donot_split_these_exons.append(start+"-"+end)
        #if len(frame_0)
        #print(len(sequence_per_exon[id])/3,len(frame_0),len(frame_1),len(frame_2))
    
    """for chromosome in donot_split_these_exons:
        print(chromosome,len(donot_split_these_exons))
        pprint.pprint(donot_split_these_exons)
    return"""
    fhr=open(exons_overlapping_with_introns_bedfilename,"r")
    for line in fhr:
        chromosome,start,end=line.strip().split()[:3]
        if chromosome not in overlaps:
            overlaps=[]
        overlaps.append(str(int(start)+1)+"-"+end)
    
    for chromosome in overlaps:
        overlaps=set(overlaps)
    
    
    #########################################################################################################
    # Read all the transcript info
    #########################################################################################################
    complete_transcript_info=readAllTranscriptsFromGTFFileInParallel([options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_transcripts_connecting_two_transcripts.gtf","dummy","dummy"])[0]
    exons_to_transcripts={}
    
    for transcript_id in complete_transcript_info:
        chromosome = complete_transcript_info[transcript_id]["chromosome"]
        if chromosome not in exons_to_transcripts:
            exons_to_transcripts[chromosome]={}
        
        for exon in complete_transcript_info[transcript_id]["exons"]:
            exon_str=str(exon[0])+"-"+str(exon[1])
            if exon_str not in exons_to_transcripts[chromosome]:
                exons_to_transcripts[chromosome][exon_str]=[]
            if transcript_id not in exons_to_transcripts[chromosome][exon_str]:
                exons_to_transcripts[chromosome][exon_str].append(transcript_id)
    
    combined_transcript_info=deepcopy(complete_transcript_info)
    
    #########################################################################################################
    # Create file for Change Points Detection computation for each condition and each exon
    #########################################################################################################
    
    # File content and structure
    # 
    # Comma separated list
    # Each entry on newline
    # Field 1 - "cov" or "opp" Indicating whether the exon should be treated for a coverage split or overlapping case
    # Field 2 - condition
    # Field 3 - Chromosome
    # Field 4 - Exon start
    # Field 5 - Exon end
    # Field 6 - Either 2 or 3 Denotes the number of break points
    # Field 7 onwards - Coverage for each nucleotide for each run
    inputfile_for_CPD=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/"+"inputfileforCPD"
    outputfile_for_CPD=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/"+"outputfileforCPD"
    
    # Create exon_coverage field
    for transcript_id in complete_transcript_info:
        complete_transcript_info[transcript_id]["exon_coverage"] = []
    
    # Populate the exon_coverage for each condition
    exons_processed_for_opp_overlapping_transcripts_all_conditions={}
    exons_processed_for_cov_transcripts_all_conditions={}
    fhw=open(inputfile_for_CPD,"w")
    for condition in options.mrna_md:
        combined_transcript_info=deepcopy(complete_transcript_info)
        for Run in options.mrna_md[condition]:
            with logging_mutex:
                logger_proxy.info("Processing inside fixOverlappingAndMergedTranscripts "+condition+" "+Run)
            coverage_info,skip=pickle.load(open(options.output_star+"/"+Run+"_counts_all_info.pkl","rb"))
            for transcript_id in coverage_info:
                chromosome=transcript_id.split(".")[0]
                temp=[]
                for exon_coverage in coverage_info[transcript_id]["bed_cov"]:
                    temp.append(np.array(exon_coverage))
                coverage_info[transcript_id]["bed_cov"]=np.array(temp)
                #if transcript_id not in complete_transcript_info:continue
                if len(complete_transcript_info[transcript_id]["exon_coverage"])==0:
                    #print("cov",coverage_info[transcript_id]["bed_cov"])
                    complete_transcript_info[transcript_id]["exon_coverage"]=np.array(coverage_info[transcript_id]["bed_cov"])
                else:
                    complete_transcript_info[transcript_id]["exon_coverage"]=np.add(complete_transcript_info[transcript_id]["exon_coverage"],
                                                                                                np.array(coverage_info[transcript_id]["bed_cov"]))
                if len(combined_transcript_info[transcript_id]["exon_coverage"])==0:
                    combined_transcript_info[transcript_id]["exon_coverage"]=np.array(coverage_info[transcript_id]["bed_cov"])
                else:
                    combined_transcript_info[transcript_id]["exon_coverage"]=np.add(combined_transcript_info[transcript_id]["exon_coverage"],
                                                                                                np.array(coverage_info[transcript_id]["bed_cov"]))
                    
        
        # Coverage patterns for each condition
        #########################################################################################################
        # Same exon shared by transcripts in opposite direction 
        #########################################################################################################
        exons_processed_for_opp_overlapping_transcripts={}
        for chromosome in exons_to_transcripts:
            exons_processed=[]
            for exon in exons_to_transcripts[chromosome]:
                if len(exons_to_transcripts[chromosome][exon])==1:continue
                list_to_be_written_to_file=[]
                for i in range(len(exons_to_transcripts[chromosome][exon])):
                    transcript_id_i=exons_to_transcripts[chromosome][exon][i]
                    j=i+1
                    while j<len(exons_to_transcripts[chromosome][exon]):
                        transcript_id_j=exons_to_transcripts[chromosome][exon][j]
                        if "+" in (complete_transcript_info[transcript_id_i]["direction"],complete_transcript_info[transcript_id_j]["direction"]) and "-" in (complete_transcript_info[transcript_id_i]["direction"],complete_transcript_info[transcript_id_j]["direction"]):
                            exon_num=0 if str(complete_transcript_info[transcript_id_i]["exons"][0][0])+"-"+str(complete_transcript_info[transcript_id_i]["exons"][0][1])==exon else -1
                            #print(exon_num,complete_transcript_info[transcript_id_i]["exon_coverage"])
                            if len(complete_transcript_info[transcript_id_i]["exon_coverage"])==0:
                                #print(condition,Run,transcript_id_i,'exon_coverage')
                                j+=1
                                continue
                            exon_coverage=complete_transcript_info[transcript_id_i]["exon_coverage"][exon_num]
                            list_to_be_written_to_file=["opp",condition,chromosome,exon.split("-")[0],exon.split("-")[1],"ext"]
                            list_to_be_written_to_file.append("2"+";"+";".join(list(map(str,exon_coverage))))
                            fhw.write(",".join(list_to_be_written_to_file)+"\n")
                            if chromosome not in exons_processed_for_opp_overlapping_transcripts:
                                exons_processed_for_opp_overlapping_transcripts=[]
                            if chromosome not in exons_processed_for_opp_overlapping_transcripts_all_conditions:
                                exons_processed_for_opp_overlapping_transcripts_all_conditions=[]
                            exons_processed_for_opp_overlapping_transcripts.append(exon)
                            exons_processed_for_opp_overlapping_transcripts_all_conditions.append(exon)
                        j+=1
                        if len(list_to_be_written_to_file)>0:
                            break
                    if len(list_to_be_written_to_file)>0:
                            break
        
        for chromosome in exons_processed_for_opp_overlapping_transcripts:
            exons_processed_for_opp_overlapping_transcripts=set(exons_processed_for_opp_overlapping_transcripts)
            
        threshold=5*len(options.mrna_md[condition])
        max_threshold=20*len(options.mrna_md[condition])
        
        exons_processed=[]
        for transcript_id in complete_transcript_info:
            # Skip mono exonic transcripts
            if len(complete_transcript_info[transcript_id]["exons"])==1:continue
            for exon_num,exon in enumerate(complete_transcript_info[transcript_id]["exons"]):
                # Skip exons that have complete ORFs in them
                # exon=str(exon[0])+"-"+str(exon[1])
                if chromosome in donot_split_these_exons and str(exon[0])+"-"+str(exon[1]) in donot_split_these_exons:
                    #print(chromosome+":"+str(exon[0])+"-"+str(exon[1]))
                    continue
                if chromosome in exons_processed_for_opp_overlapping_transcripts and str(exon[0])+"-"+str(exon[1]) in exons_processed_for_opp_overlapping_transcripts:continue    
                if str(exon[0])+"-"+str(exon[1]) not in exons_processed:
                    exons_processed.append(str(exon[0])+"-"+str(exon[1]))
                else:
                    continue
                # Skip exons that overlap with introns since coverage will dip anyways 
                if chromosome in overlaps and str(exon[0])+"-"+str(exon[1]) in overlaps:continue
                external_exon=(1 if (exon_num==0 or exon_num==len(complete_transcript_info[transcript_id]["exons"])-1) else 0)
                try:
                    exon_coverage=complete_transcript_info[transcript_id]["exon_coverage"][exon_num]
                except IndexError:
                    print(transcript_id,len(complete_transcript_info[transcript_id]["exons"]),len(complete_transcript_info[transcript_id]["exon_coverage"]))
                    sys.exit()
                if external_exon==0:
                    min_coverage=min(exon_coverage)
                    max_coverage=max(exon_coverage)
                else:
                    # Skip the external exon if its less than 100
                    if len(exon_coverage)<100:continue
                    if exon_num==0: # First Exon
                        #start_from=100 if len(exon_coverage)>100 else int(0.1*len(exon_coverage))
                        min_coverage=min(exon_coverage)
                        max_coverage=max(exon_coverage)
                    else: # Last Exon
                        #end_at=100 if len(exon_coverage)>100 else int(0.1*len(exon_coverage))
                        min_coverage=min(exon_coverage)
                        max_coverage=max(exon_coverage)
                std_dev=np.std(exon_coverage)
                
                if min_coverage<threshold and std_dev>=5 and len(exon_coverage)>=50:
                    list_to_be_written_to_file=["cov",condition,str(chromosome),str(exon[0]),str(exon[1]),("ext" if external_exon==1 else "int")]
                    #list_to_be_written_to_file.extend(list(map(str,exon_coverage)))
                    list_to_be_written_to_file.append(("3" if external_exon==1 else "2")+";"+";".join(list(map(str,exon_coverage))))
                    #all_lines_to_be_written["cov"][str(chromosome)+":"+str(exon[0])+"-"+str(exon[1])]=list_to_be_written_to_file
                    fhw.write(",".join(list_to_be_written_to_file)+"\n")
                    if chromosome not in exons_processed_for_cov_transcripts_all_conditions:
                        exons_processed_for_cov_transcripts_all_conditions=[]
                    exons_processed_for_cov_transcripts_all_conditions.append(str(exon[0])+"-"+str(exon[1]))
        # Clearing out complete_transcript_info for next run
        for transcript_id in complete_transcript_info:
            complete_transcript_info[transcript_id]["exon_coverage"]=[]
    
    for chromosome in exons_processed_for_opp_overlapping_transcripts_all_conditions:
        exons_processed_for_opp_overlapping_transcripts_all_conditions=set(exons_processed_for_opp_overlapping_transcripts_all_conditions)
    
    for chromosome in exons_processed_for_cov_transcripts_all_conditions:
        exons_processed_for_cov_transcripts_all_conditions=set(exons_processed_for_cov_transcripts_all_conditions)
    
    # Coverage patterns for all conditions combined
    exons_in_opp=[]
    exons_in_cov=[]
    
    with logging_mutex:
        logger_proxy.info("Processing inside fixOverlappingAndMergedTranscripts all conditions and runs combined")
    for transcript_id in combined_transcript_info:
        chromosome=combined_transcript_info[transcript_id]["chromosome"]
        for exon_num,exon in enumerate(combined_transcript_info[transcript_id]["exons"]):
            exon_str=str(exon[0])+"-"+str(exon[1])
            if chromosome in exons_processed_for_opp_overlapping_transcripts_all_conditions and exon_str in exons_processed_for_opp_overlapping_transcripts_all_conditions[chromosome]:
                if exon_str not in exons_in_opp:
                    exons_in_opp.append(exon_str)
                else:
                    continue
                exon_coverage=combined_transcript_info[transcript_id]["exon_coverage"][exon_num]
                list_to_be_written_to_file=["opp","combined",chromosome,exon_str.split("-")[0],exon_str.split("-")[1],"ext"]
                list_to_be_written_to_file.append("2"+";"+";".join(list(map(str,exon_coverage))))
                fhw.write(",".join(list_to_be_written_to_file)+"\n")
                
            if chromosome in exons_processed_for_cov_transcripts_all_conditions and exon_str in exons_processed_for_cov_transcripts_all_conditions:
                if exon_str not in exons_in_cov:
                    exons_in_cov.append(exon_str)
                else:
                    continue
                if exon_str in exons_in_opp:
                    print(chromosome,exon,"found in opp")
                    sys.stdout.flush()
                    continue
                external_exon=(1 if (exon_num==0 or exon_num==len(combined_transcript_info[transcript_id]["exons"])-1) else 0)
                exon_coverage=combined_transcript_info[transcript_id]["exon_coverage"][exon_num]
                list_to_be_written_to_file=["cov","combined",str(chromosome),str(exon[0]),str(exon[1]),("ext" if external_exon==1 else "int")]
                #list_to_be_written_to_file.extend(list(map(str,exon_coverage)))
                list_to_be_written_to_file.append(("3" if external_exon==1 else "2")+";"+";".join(list(map(str,exon_coverage))))
                #all_lines_to_be_written["cov"][str(chromosome)+":"+str(exon[0])+"-"+str(exon[1])]=list_to_be_written_to_file
                fhw.write(",".join(list_to_be_written_to_file)+"\n")
    
    inputfile_for_CPD=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/"+"inputfileforCPD"
    outputfile_for_CPD=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/"+"outputfileforCPD"
    pool = multiprocessing.Pool(processes=int(options.cpu))
    # Split inputfile into chunks
    if len(open(inputfile_for_CPD,"r").read().split("\n"))>5000:
        chunked_file=list(divide_chunks(open(inputfile_for_CPD,"r").read().split("\n"),5000))
        time.sleep(5)
        
        for file_num,chunk in enumerate(chunked_file):
            fhw=open(inputfile_for_CPD+"_"+str(file_num),"w")
            fhw.write("\n".join(chunk))
            fhw.close()
    else:
        os.system(f"cp {inputfile_for_CPD} {inputfile_for_CPD}_0")
    
    with logging_mutex:
        logger_proxy.info("inputfileforCPD generation is complete")
    
    round_num=1
    present=0
    while present==0:
        with logging_mutex:
            logger_proxy.info("Starting round "+str(round_num))
        allinputs=[]
        for file_num,chunk in enumerate(chunked_file):
            outputfilename_portion=outputfile_for_CPD+"_"+str(file_num)
            cmd="Rscript "+options.softwares["find_exonic_troughs"]
            cmd+=" "+inputfile_for_CPD+"_"+str(file_num)
            cmd+=" "+outputfilename_portion
            cmd+=" 2 "
            cmd+=" 2> "+outputfilename_portion+".error"
            
            if os.path.exists(outputfilename_portion)==False:
                allinputs.append(["dummy",cmd]) 
                with logging_mutex:
                    logger_proxy.info("Performing calculation for "+outputfile_for_CPD+"_"+str(file_num))
        
        pool.map(runCommand,allinputs)
        time.sleep(5)
         
        file_absent=0
        for file_num,chunk in enumerate(chunked_file):
            if os.path.exists(outputfile_for_CPD+"_"+str(file_num))==True:
                pass
            else:
                file_absent=1
                with logging_mutex:
                    logger_proxy.info("Missing outputfile "+outputfile_for_CPD+"_"+str(file_num))
        if file_absent==1:
            present=0
        elif file_absent==0:
            present=1
        if present==0: 
            round_num+=1
    
    # Combine results in one file
    cmd="cat "
    for f in range(len(list(chunked_file))):
        cmd+=outputfile_for_CPD+"_"+str(f)+" "
    cmd+=" > "+outputfile_for_CPD
    cmd+=" 2> "+outputfile_for_CPD+".error"
    os.system(cmd)
    #return
    cmd="rm "+outputfile_for_CPD+"_* "+inputfile_for_CPD+"_* "
    os.system(cmd)
    
    inputfile_for_CPD=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/"+"inputfileforCPD"
    outputfile_for_CPD=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/"+"outputfileforCPD"
    
    fhr_coverage=open(inputfile_for_CPD,"r")
    fhr_changepoints=open(outputfile_for_CPD,"r")
    exons_with_no_breakpoints_opp={}
    exons_to_split_exons_opp={}
    exons_to_split_exons_cov={}
    complete_info_about_breakpoints={"cov":{},"opp":{}}    
    while True:
        line_coverage=fhr_coverage.readline().strip()
        line_changepoints=fhr_changepoints.readline().strip()
        if not line_coverage or not line_changepoints:
            break
        type_cov,condition_cov,chromosome_cov,start_cov,end_cov,exon_loc_cov,coverage=line_coverage.strip().split(",")
        type_cpd,condition_cpd,chromosome_cpd,start_cpd,end_cpd,num_breakpoints,exon_loc_cpd,changepoints_cpd=line_changepoints.strip().split(",")
        num_breakpoints=int(num_breakpoints)
        if ";" in changepoints_cpd:
            changepoints_cpd=list(map(int,changepoints_cpd.split(";")))
        else:
            changepoints_cpd=[]
        coverage=list(map(int,coverage.split(";")))[2:]
        exon=start_cpd+"-"+end_cpd
        
        if type_cov=="opp":
            if len(changepoints_cpd)==0:
                pass
            else:
                if chromosome_cpd not in exons_to_split_exons_opp:
                    exons_to_split_exons_opp[chromosome_cpd]={}
                if len(changepoints_cpd)==1:
                    exons_to_split_exons_opp[chromosome_cpd][start_cpd+"-"+end_cpd]=[int(changepoints_cpd[0])]
                elif len(changepoints_cpd)==2:
                    exons_to_split_exons_opp[chromosome_cpd][start_cpd+"-"+end_cpd]=[int(changepoints_cpd[0]),int(changepoints_cpd[1])]
                if chromosome_cpd not in complete_info_about_breakpoints["opp"]:
                    complete_info_about_breakpoints["opp"][chromosome_cpd]={}
                if exon not in complete_info_about_breakpoints["opp"][chromosome_cpd]:
                    complete_info_about_breakpoints["opp"][chromosome_cpd][exon]={}
                complete_info_about_breakpoints["opp"][chromosome_cpd][exon][condition_cpd]={"changepoints":changepoints_cpd,"ratio":[]}
        elif type_cov=="cov": 
            if chromosome_cpd not in exons_to_split_exons_cov:
                exons_to_split_exons_cov[chromosome_cpd]={}
            if num_breakpoints==1:
                pass
            elif num_breakpoints==2:
                if len(changepoints_cpd)==2:
                    coverage=list(map(int,coverage))
                    coverage_1st_portion=coverage[:changepoints_cpd[0]]
                    coverage_2nd_portion=coverage[changepoints_cpd[0]:changepoints_cpd[1]]
                    coverage_3rd_portion=coverage[changepoints_cpd[1]:]
                    ratio1=round(np.average(coverage_2nd_portion)/np.average(coverage_1st_portion),2)
                    ratio2=round(np.average(coverage_2nd_portion)/np.average(coverage_3rd_portion),2)
                    """
                    if np.average(coverage_2nd_portion)/np.average(coverage_1st_portion)<1 and np.average(coverage_2nd_portion)/np.average(coverage_3rd_portion)<1:
                        print(condition,chromosome_cov+":"+start_cov+"-"+end_cov,np.average(coverage_2nd_portion)/np.average(coverage_1st_portion),np.average(coverage_2nd_portion)/np.average(coverage_3rd_portion),changepoints)
                    """
                    exons_to_split_exons_cov[chromosome_cpd][start_cpd+"-"+end_cpd]=[int(changepoints_cpd[0]),int(changepoints_cpd[1])]
                    #print(line_changepoints,ratio1,ratio2)
                    if 0<min(ratio1,ratio2)<1 and max(ratio1,ratio2)<1:
                        if chromosome_cpd not in complete_info_about_breakpoints["cov"]:
                            complete_info_about_breakpoints["cov"][chromosome_cpd]={}
                        if exon not in complete_info_about_breakpoints["cov"][chromosome_cpd]:
                            complete_info_about_breakpoints["cov"][chromosome_cpd][exon]={}
                        complete_info_about_breakpoints["cov"][chromosome_cpd][exon][condition_cpd]={"changepoints":changepoints_cpd,"ratio":[ratio1,ratio2]}
                    
                    """if chromosome_cpd=="4" and start_cpd=="10439766" and end_cpd=="10440165":
                        print("Ratios",ratio1,ratio2)
                        sys.stdout.flush()"""
                else:
                    # Do nothing - its a false positive
                    if chromosome_cpd not in complete_info_about_breakpoints["cov"]:
                        complete_info_about_breakpoints["cov"][chromosome_cpd]={}
                    if exon not in complete_info_about_breakpoints["cov"][chromosome_cpd]:
                        complete_info_about_breakpoints["cov"][chromosome_cpd][exon]={}
                    complete_info_about_breakpoints["cov"][chromosome_cpd][exon][condition_cpd]={"changepoints":[],"ratio":[]}
                
            elif num_breakpoints==3:
                if len(changepoints_cpd)==3:
                    coverage=list(map(int,coverage))
                    coverage_1st_portion=coverage[:changepoints_cpd[0]]
                    coverage_2nd_portion=coverage[changepoints_cpd[0]:changepoints_cpd[1]]
                    coverage_3rd_portion=coverage[changepoints_cpd[1]:changepoints_cpd[2]]
                    coverage_4th_portion=coverage[changepoints_cpd[2]:]
                    
                    ratio1=round(np.average(coverage_2nd_portion)/np.average(coverage_1st_portion),2)
                    ratio2=round(np.average(coverage_2nd_portion)/np.average(coverage_3rd_portion),2)
                    ratio3=round(np.average(coverage_3rd_portion)/np.average(coverage_2nd_portion),2)
                    ratio4=round(np.average(coverage_3rd_portion)/np.average(coverage_4th_portion),2)
                    
                    if 0<min(ratio1,ratio2)<1 and max(ratio1,ratio2)<1:
                        #print(condition,changepoints,chromosome_cov+":"+start_cov+"-"+end_cov,np.average(coverage_2nd_portion)/np.average(coverage_1st_portion),np.average(coverage_2nd_portion)/np.average(coverage_3rd_portion))
                        #print(line_changepoints,ratio1,ratio2)
                        if chromosome_cpd not in complete_info_about_breakpoints["cov"]:
                            complete_info_about_breakpoints["cov"][chromosome_cpd]={}
                        if exon not in complete_info_about_breakpoints["cov"][chromosome_cpd]:
                            complete_info_about_breakpoints["cov"][chromosome_cpd][exon]={}
                        complete_info_about_breakpoints["cov"][chromosome_cpd][exon][condition_cpd]={"changepoints":changepoints_cpd[:2],"ratio":[ratio1,ratio2]}
                    elif 0<min(ratio3,ratio4)<1 and max(ratio3,ratio4)<1:
                        #print(condition,changepoints,chromosome_cov+":"+start_cov+"-"+end_cov,np.average(coverage_3rd_portion)/np.average(coverage_2nd_portion),np.average(coverage_3rd_portion)/np.average(coverage_4th_portion))
                        #print(line_changepoints,ratio3,ratio4)
                        if chromosome_cpd not in complete_info_about_breakpoints["cov"]:
                            complete_info_about_breakpoints["cov"][chromosome_cpd]={}
                        if exon not in complete_info_about_breakpoints["cov"][chromosome_cpd]:
                            complete_info_about_breakpoints["cov"][chromosome_cpd][exon]={}
                        complete_info_about_breakpoints["cov"][chromosome_cpd][exon][condition_cpd]={"changepoints":changepoints_cpd[1:],"ratio":[ratio3,ratio4]}
                    """if chromosome_cpd=="5" and start_cpd=="6688611" and end_cpd=="6690038":
                        print("this",ratio1,ratio2,ratio3,ratio4,np.average(coverage_1st_portion),np.average(coverage_2nd_portion),np.average(coverage_3rd_portion),np.average(coverage_4th_portion))"""
                else:
                    # Do nothing - its a false positive
                    if chromosome_cpd not in complete_info_about_breakpoints["cov"]:
                        complete_info_about_breakpoints["cov"][chromosome_cpd]={}
                    if exon not in complete_info_about_breakpoints["cov"][chromosome_cpd]:
                        complete_info_about_breakpoints["cov"][chromosome_cpd][exon]={}
                    complete_info_about_breakpoints["cov"][chromosome_cpd][exon][condition_cpd]={"changepoints":[],"ratio":[]}
    
    # Modify end exons of transcripts
    overlapping_transcripts_end_modified=[]
    for chromosome in complete_info_about_breakpoints["opp"]:        
        for exon in complete_info_about_breakpoints["opp"][chromosome]:
            #if len(complete_info_about_breakpoints["cov"][chromosome][exon])==total_num_of_conditions+1:
            median_changepoints=list(map(int,list(np.median(np.array([complete_info_about_breakpoints["opp"][chromosome][exon][condition]["changepoints"] for condition in complete_info_about_breakpoints["opp"][chromosome][exon]]),axis=0))))
            for transcript_id in exons_to_transcripts[chromosome][exon]:
                #print("Altering",transcript_id)
                #print(exon)
                exon_list=list(map(int,exon.split("-")))
                #pprint.pprint(complete_transcript_info[transcript_id]["exons"])
                if complete_transcript_info[transcript_id]["exons"][0][0]==exon_list[0] and complete_transcript_info[transcript_id]["exons"][0][1]==exon_list[1]:
                    #print("Inside this 1")
                    complete_transcript_info[transcript_id]["exons"][0][0]=complete_transcript_info[transcript_id]["exons"][0][0]+median_changepoints[0]
                elif complete_transcript_info[transcript_id]["exons"][-1][0]==exon_list[0] and complete_transcript_info[transcript_id]["exons"][-1][1]==exon_list[1]:
                    #print("Inside this 2")
                    complete_transcript_info[transcript_id]["exons"][-1][1]=complete_transcript_info[transcript_id]["exons"][-1][0]+median_changepoints[1]
                overlapping_transcripts_end_modified.append(transcript_id)
                """pprint.pprint(complete_transcript_info[transcript_id]["exons"])
                sys.stdout.flush()"""
                
    # Select changepoints for coverage based on the conditions they were recognized in
    total_num_of_conditions=len(options.mrna_md)
    #exon_to_cov_breakpoints={}
    all_cov_breakpoints=[]
    for chromosome in complete_info_about_breakpoints["cov"]:
        #exon_to_cov_breakpoints[chromosome]={}
        for exon in complete_info_about_breakpoints["cov"][chromosome]:
            #if len(complete_info_about_breakpoints["cov"][chromosome][exon])==total_num_of_conditions+1:
            collect_changepoints=[]
            median_changepoints=[]
            tried_in_num_conditions=len([1 for condition in complete_info_about_breakpoints["cov"][chromosome][exon]])
            
            for condition in complete_info_about_breakpoints["cov"][chromosome][exon]:
                if condition=="combined":continue
                if len(complete_info_about_breakpoints["cov"][chromosome][exon][condition]["changepoints"])>0:
                    try:
                        collect_changepoints.append(complete_info_about_breakpoints["cov"][chromosome][exon][condition]["changepoints"])
                        median_changepoints=list(np.median(np.array(collect_changepoints),axis=0))
                    except TypeError:
                        print("ERROR")
                        print(collect_changepoints)
                        sys.exit()
            if len(collect_changepoints)==0:
                continue
            if len(collect_changepoints)/tried_in_num_conditions<0.50:continue
            if len(median_changepoints)==0:continue
            #median_changepoints=list(np.median(np.array([complete_info_about_breakpoints["cov"][chromosome][exon][condition]["changepoints"] for condition in complete_info_about_breakpoints["cov"][chromosome][exon]]),axis=0))
            """
            print(chromosome+":"+exon,len(complete_info_about_breakpoints["cov"][chromosome][exon]),"cov",[condition for condition in complete_info_about_breakpoints["cov"][chromosome][exon]],
                  [(complete_info_about_breakpoints["cov"][chromosome][exon][condition]["changepoints"],complete_info_about_breakpoints["cov"][chromosome][exon][condition]["ratio"]) for condition in complete_info_about_breakpoints["cov"][chromosome][exon]],
                  median_changepoints)
            """
            #if len(complete_info_about_breakpoints["cov"][chromosome][exon])==1:continue
            exon=list(map(int,exon.split("-")))
            #exon_to_cov_breakpoints[chromosome][exon]=median_changepoints
            try:
                all_cov_breakpoints.append([chromosome,exon[0]+int(median_changepoints[0]),exon[0]+int(median_changepoints[1])])
            except IndexError:
                print("ERROR")
                print(chromosome,exon[0],exon[1],median_changepoints,collect_changepoints)
                sys.exit()
            #print(chromosome,exon[0]+int(median_changepoints[0]),"-",exon[0]+int(median_changepoints[1]))
    
    #pprint.pprint(complete_info_about_breakpoints["cov"])
    all_cov_breakpoints_pd=pd.DataFrame(all_cov_breakpoints,columns=["chromosome","cutpoint1","cutpoint2"])
    all_cov_breakpoints_pd = all_cov_breakpoints_pd.sort_values(["chromosome","cutpoint1"])
    """with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(all_cov_breakpoints_pd)"""
    
    all_transcript_info=[]
    for transcript_id in complete_transcript_info:
        exon_definition=complete_transcript_info[transcript_id]["exons"]
        transcript_start=complete_transcript_info[transcript_id]["exons"][0][0]
        transcript_end=complete_transcript_info[transcript_id]["exons"][-1][1]
        direction=complete_transcript_info[transcript_id]["direction"]
        all_transcript_info.append([chromosome,transcript_id,exon_definition,transcript_start,transcript_end,direction])
    all_transcript_info_pd=pd.DataFrame(all_transcript_info,columns=["chromosome","transcript_id","exon_definition","transcript_start","transcript_end","direction"])
    
    prev_to_new_transcripts={}
    remove_these_transcripts=[]
    exons_affected_by_coverage_splits={}
    for index,row in all_transcript_info_pd.iterrows():
        chromosome,transcript_id,exon_definition,transcript_start,transcript_end,direction=row
        selected_breakpoints=all_cov_breakpoints_pd[(all_cov_breakpoints_pd["chromosome"]==chromosome) & (transcript_start<=all_cov_breakpoints_pd["cutpoint1"]) & (all_cov_breakpoints_pd["cutpoint2"]<=transcript_end)]
        if selected_breakpoints.shape[0]==0:continue
        #print("transcript_id",transcript_id)
        breakpoints_list_of_lists=[]
        for index,breakpoint_row in selected_breakpoints.iterrows():
            chromosome,cutpoint1,cutpoint2=breakpoint_row
            breakpoints_list_of_lists.append([cutpoint1,cutpoint2,'c'])
        
        exon_list_of_lists=[]
        for exon in exon_definition:
            exon_list_of_lists.append([exon[0],exon[1],'e'])
        
        combined_list_of_exons_and_breakpoints=[]
        combined_list_of_exons_and_breakpoints.extend(breakpoints_list_of_lists)
        combined_list_of_exons_and_breakpoints.extend(exon_list_of_lists)
        #pprint.pprint(combined_list_of_exons_and_breakpoints)
        #print("-"*100)
        combined_list_of_exons_and_breakpoints.sort(key = lambda x: x[0])
        #pprint.pprint(combined_list_of_exons_and_breakpoints)
        #print("="*100) 
        if chromosome not in exons_affected_by_coverage_splits:
            exons_affected_by_coverage_splits={}
            
        exon_definition_of_split_transcripts=createNewTranscripts(combined_list_of_exons_and_breakpoints,transcript_id)
        remove_these_transcripts.append(transcript_id)
        if transcript_id not in prev_to_new_transcripts:
            prev_to_new_transcripts[transcript_id]=[]
            
        for new_transcript_num,each_new_transcript_exon_definition in enumerate(exon_definition_of_split_transcripts):
            gene_id=".".join(transcript_id.split(".")[:2])
            transcript_num=transcript_id.split(".")[-1]
            new_transcript_id=gene_id+"_"+transcript_num+"_covsplit"+"."+str(new_transcript_num)
            chromosome=new_transcript_id.split(".")[0]
            cov=complete_transcript_info[transcript_id]["cov"]
            tpm=complete_transcript_info[transcript_id]["TPM"]
            fpkm=complete_transcript_info[transcript_id]["FPKM"]
            complete_transcript_info[new_transcript_id]={"exons":[],
                                                                 "introns":[],
                                                                 "exon_coverage":[],
                                                                 "direction":direction  ,
                                                                 "exon_number_overlapping_with_introns_0_based":[],
                                                                 "cov":cov,
                                                                 "TPM":tpm,
                                                                 "FPKM":fpkm,
                                                                 "chromosome":chromosome,
                                                                 "annotator":"FINDER"
                                                                 }
            complete_transcript_info[new_transcript_id]["exons"]=each_new_transcript_exon_definition
            prev_to_new_transcripts[transcript_id].append(new_transcript_id)
            i=1
            while i<len(complete_transcript_info[new_transcript_id]["exons"]):
                complete_transcript_info[new_transcript_id]["introns"].append([complete_transcript_info[new_transcript_id]["exons"][i-1][1]+1,complete_transcript_info[new_transcript_id]["exons"][i][0]-1])
                i+=1
                
    # Remove the presplit transcripts
    for transcript_id in remove_these_transcripts:
        del complete_transcript_info[transcript_id]
    
    complete_transcript_info = removeOverlappingExonsFromEachTranscript(complete_transcript_info)
    #complete_transcript_info = removeSpuriousExonsAndTranscripts(complete_transcript_info)
    
    writeTranscriptsToFile([complete_transcript_info,output_gtf_filename,0])
    cmd=f"cat {output_gtf_filename}|uniq > {output_gtf_filename}.temp"
    os.system(cmd)
    cmd=f"mv {output_gtf_filename}.temp {output_gtf_filename}"
    os.system(cmd)
    
    