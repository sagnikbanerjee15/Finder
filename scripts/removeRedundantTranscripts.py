#! /usr/bin/env python

from scripts.fileReadWriteOperations import *
import multiprocessing
import os
import sys

import pandas as pd


def findSubsetTranscripts(eachinput):
    """
    """
    gene_to_transcripts,gene,transcripts_fasta=eachinput
    redundant_transcripts=[]
    i=0
    while i<len(gene_to_transcripts[gene]):
        transcript_i=gene_to_transcripts[gene][i]
        j=i+1
        while j<len(gene_to_transcripts[gene]):
            transcript_j=gene_to_transcripts[gene][j]
            if transcripts_fasta[transcript_i] in transcripts_fasta[transcript_j]:
                redundant_transcripts.append(transcript_i)
            if transcripts_fasta[transcript_j] in transcripts_fasta[transcript_i]:
                redundant_transcripts.append(transcript_j)
            j+=1
        i+=1
    
    return gene,list(set(redundant_transcripts))
    
def removeRedundantTranscripts(input_gtf_filename,output_gtf_filename,options):
    if os.path.exists(output_gtf_filename)==True:return
    redundant_transcripts=[]
    gene_to_transcripts={}
    transcript_all_info=readAllTranscriptsFromGTFFileInParallel([input_gtf_filename,"dummy","dummy"])[0]
    
    # Remove exact same copies of transcripts from whole genome
    
    
    remove_these=[]
    exon_definitons=[]
    for transcript_id in transcript_all_info:
        all_exons="_".join([str(exon[0])+"-"+str(exon[1]) for exon in transcript_all_info[transcript_id]["exons"]])
        if all_exons in set(exon_definitons):
            remove_these.append(transcript_id)
            redundant_transcripts.append(transcript_id)
        else:
            exon_definitons.append(all_exons)
    #print("\n".join(remove_these))
    for transcript_id in remove_these:
        del transcript_all_info[transcript_id]
    
    writeTranscriptsToFile([transcript_all_info,input_gtf_filename+".temp",0])
    transcript_all_info=readAllTranscriptsFromGTFFileInParallel([input_gtf_filename+".temp","dummy","dummy"])[0]
    
    for transcript_id in transcript_all_info:
        gene=".".join(transcript_id.split(".")[:2])
        if gene not in gene_to_transcripts:
            gene_to_transcripts[gene]=[]
        gene_to_transcripts[gene].append(transcript_id)
    
    cmd="gffread "
    cmd+=" -w "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined.fasta "
    cmd+=" -g "+options.genome+" "
    cmd+=input_gtf_filename+".temp "
    os.system(cmd)
    
    cmd="perl -pe '/^>/ ? print \"\n\" : chomp' "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined.fasta "
    cmd+="| tail -n +2 > "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_sl.fasta"
    os.system(cmd)
    
    cmd="mv "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_sl.fasta "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined.fasta"
    os.system(cmd)
    
    transcripts_fasta=readFastaFile(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined.fasta")
    pool = multiprocessing.Pool(processes=int(options.cpu))
    all_inputs=[]
    for gene in gene_to_transcripts:
        all_inputs.append([gene_to_transcripts,gene,transcripts_fasta])
    
    results=pool.map(findSubsetTranscripts,all_inputs)
    
    for result in results:
        redundant_transcripts.extend(result[-1])
    
    #print("\n".join(redundant_transcripts))
    # Remove uniexon transcripts
    #complete_transcript_info,useless1,useless2=readAllTranscriptsFromGTFFileInParallel([input_gtf_filename,"dummy","dummy"])
    complete_transcript_info=transcript_all_info
    
    all_transcript_info=[]
    uni_exon_transcript_info=[]
    
    for transcript_id in complete_transcript_info:
        chromosome = complete_transcript_info[transcript_id]["chromosome"]
        exon_definition=complete_transcript_info[transcript_id]["exons"]
        transcript_start=complete_transcript_info[transcript_id]["exons"][0][0]
        transcript_end=complete_transcript_info[transcript_id]["exons"][-1][1]
        direction=complete_transcript_info[transcript_id]["direction"]
        if len(exon_definition)>1:
            for exon in exon_definition:
                all_transcript_info.append([chromosome,transcript_id,exon[0],exon[1],transcript_start,transcript_end,direction])
        else:
            for exon in exon_definition:
                uni_exon_transcript_info.append([chromosome,transcript_id,exon[0],exon[1],transcript_start,transcript_end,direction])
                
    all_transcript_info_pd=pd.DataFrame(all_transcript_info,columns=["chromosome","transcript_id","exon_start","exon_end","transcript_start","transcript_end","direction"])
    uni_exon_transcripts_pd=pd.DataFrame(uni_exon_transcript_info,columns=["chromosome","transcript_id","exon_start","exon_end","transcript_start","transcript_end","direction"])
    for row_num,row in uni_exon_transcripts_pd.iterrows():
        chromosome,uni_exon_transcript,exon_start,exon_end,transcript_start,transcript_end,direction=row
        overlapping_exons_pd=all_transcript_info_pd[(all_transcript_info_pd["chromosome"]==chromosome) & (all_transcript_info_pd["exon_start"]<=exon_start) & (all_transcript_info_pd["exon_end"]>=exon_end)]
        if overlapping_exons_pd.shape[0]>0:
            for row_num,row in overlapping_exons_pd.iterrows():
                chromosome,each_redundant_transcript,exon_start,exon_end,transcript_start,transcript_end,direction=row
                redundant_transcripts.append(uni_exon_transcript)
    
    # Find uni-exon transcripts which are exactly same or subsets of other uni-exon transcripts
    uni_exon_transcripts_pd.sort_values(by=['chromosome', 'transcript_start'])
    uni_exon_transcript_info=uni_exon_transcripts_pd.values.tolist()
    i=0
    while i<len(uni_exon_transcript_info):
        j=i+1
        chromosome_i,uni_exon_transcript_i,exon_start_i,exon_end_i,transcript_start_i,transcript_end_i,direction_i=uni_exon_transcript_info[i]
        if (transcript_end_i-transcript_start_i)<150:
            if uni_exon_transcript_i not in set(redundant_transcripts):
                redundant_transcripts.append(uni_exon_transcript_i)
                #print("this","Length of redundant transcripts",len(redundant_transcripts),"Number of redundant transcripts",len(set(redundant_transcripts)))
                sys.stdout.flush()
            i+=1
            continue
        while j<len(uni_exon_transcript_info):
            chromosome_j,uni_exon_transcript_j,exon_start_j,exon_end_j,transcript_start_j,transcript_end_j,direction_j=uni_exon_transcript_info[j]
            if chromosome_j!=chromosome_i:break
            if transcript_start_j>transcript_end_i:break
            if transcript_start_i<=transcript_start_j<=transcript_end_j<=transcript_end_i or transcript_start_j<=transcript_start_i<=transcript_end_i<=transcript_end_j:
                #print("this",i,j,uni_exon_transcript_i,uni_exon_transcript_j,transcript_start_i,transcript_end_i,transcript_start_j,transcript_end_j)
                if uni_exon_transcript_i not in set(redundant_transcripts):
                    redundant_transcripts.append(uni_exon_transcript_i)
                    #print("this","Length of redundant transcripts",len(redundant_transcripts),"Number of redundant transcripts",len(set(redundant_transcripts)))
                    sys.stdout.flush()
                break
            j+=1
        i+=1
    

    redundant_transcripts=set(redundant_transcripts)
    #print("\n".join(redundant_transcripts))
    sys.stdout.flush()
    fhr=open(input_gtf_filename,"r")
    fhw=open(output_gtf_filename,"w")
    for line in fhr:
        if line[0]=="#":continue
        chromosome,psiclass,structure,start,end,useless1,direction,useless2,desc=line.strip().split("\t")
        exon=start+"-"+end
        for ele in desc.split(";"):
            if "gene_id" in ele:
                gene_id=ele.split()[-1].strip("\"")
            if "transcript_id" in ele:
                transcript_id=ele.split()[-1].strip("\"")
        if transcript_id in redundant_transcripts:
            continue
        else:
            fhw.write(line)
    fhr.close()
    fhw.close()
    os.system("rm -rf "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined.fasta")
