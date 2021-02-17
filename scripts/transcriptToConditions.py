

from scripts.fileReadWriteOperations import *
import collections
import multiprocessing


def transcriptToConditions(options):
    """
    """
    output_filename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/transcript_to_condition"
    transcript_to_condition={}
    pool = multiprocessing.Pool(processes=int(options.cpu))
    #filename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_merged_transcripts.gtf"
    filename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_split_transcripts_with_bad_SJ_redundancy_removed.gtf"
    final_transcripts_temp,useless1,useless2=readAllTranscriptsFromGTFFileInParallel([filename,"dummy","dummy"])
    final_transcripts=[transcript_id for transcript_id in final_transcripts_temp]
    
    transcripts_from_individual_runs={}
    allinputs=[]
    for condition in options.mrna_md:
        transcripts_from_individual_runs[condition]=[]
        for Run in options.mrna_md[condition]:
            filename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/psiclass_output_sample_"+Run+".gtf"
            allinputs.append([filename,Run,condition])
    
    results=pool.map(readAllTranscriptsFromGTFFileInParallel,allinputs)
    for result in results:
        info,Run,condition=result
        transcripts_from_individual_runs[condition].extend([transcript_id for transcript_id in info])
        
    for condition in transcripts_from_individual_runs:
        transcripts_from_individual_runs[condition]=dict(collections.Counter(transcripts_from_individual_runs[condition]))
    
    for transcript_id in final_transcripts:
        transcript_to_condition[transcript_id]=[]
        for condition in transcripts_from_individual_runs:
            if "NODE" in transcript_id:continue
            if "cov" in transcript_id:
                transcript_id_new=transcript_id.split("_cov")[0]
                transcript_id_new=transcript_id_new.replace("_",".")
                #print(transcript_id,transcript_id_new,transcript_id_new in transcripts_from_individual_runs[condition],condition)
                if transcript_id_new in transcripts_from_individual_runs[condition]:
                    transcript_to_condition[transcript_id].append(condition)
            elif "merged" in transcript_id:
                transcript_id_new1,transcript_id_new2=transcript_id.split("_merged_")
                #print(transcript_id,transcript_id_new1,transcript_id_new2,transcript_id_new1 in transcripts_from_individual_runs[condition] or transcript_id_new2 in transcripts_from_individual_runs[condition],condition)
                if transcript_id_new1 in transcripts_from_individual_runs[condition] or transcript_id_new2 in transcripts_from_individual_runs[condition]==True:
                    transcript_to_condition[transcript_id].append(condition)
            else:
                #print(transcript_id,transcript_id in transcripts_from_individual_runs[condition],condition)
                if transcript_id in transcripts_from_individual_runs[condition]:
                    transcript_to_condition[transcript_id].append(condition)
    fhw=open(output_filename,"w")
    for transcript_id in transcript_to_condition:
        #print(transcript_id,transcript_to_condition[transcript_id])
        fhw.write(transcript_id+"\t"+",".join(transcript_to_condition[transcript_id])+"\n")
    fhw.close()
   
