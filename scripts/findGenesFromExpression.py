##############################################################################################################################################################################
# 
# Program contains functions to construct gene structures from expressed RNA-Seq data 
# 
##############################################################################################################################################################################

from scripts.alignReads import *
from scripts.fileReadWriteOperations import *
import multiprocessing
import os
import sys
import random


def selectHighConfidenceSpliceJunctionsPerCondition(eachinput):

    condition,options,round=eachinput
    
    # Merge SJ.out.tab files from all samples of the same condition
    
    sjdbout_filename_merged_condition=options.output_star+"/"+condition+"_round"+str(round)+"_SJ.out.tab"
    
    num_samples=0
    cmd="cat "
    for Run in options.mrna_md[condition]:
        sjdbout_filename_round1=options.output_star+"/"+Run+"_round"+str(round)+"_SJ.out.tab "
        cmd+=sjdbout_filename_round1
        num_samples+=1
    if round==2:
        cmd+=" |awk '$6==0' |sort > "+sjdbout_filename_merged_condition
    else:
        cmd+=" |sort > "+sjdbout_filename_merged_condition
    os.system(cmd)
    #print(cmd)
    
    junctions_to_be_discarded={}
    junctions_to_be_retained={}
    fhr=open(sjdbout_filename_merged_condition,"r")
    prev_chromosome,prev_junction="",""
    same_junction=[]
    #fraction_of_samples_junction_should_be_in=1 if num_samples==1 else (1/2 if num_samples==2 else 2/3)
    prev_line=""
    while True:
        line=fhr.readline()
        if not line:break
        chromosome,j_start,j_end,strand,intron_motif,annotated,uniq_reads,mm_reads,max_overhang=line.strip().split()
        if chromosome not in junctions_to_be_discarded:
            junctions_to_be_discarded[chromosome]=[]
            junctions_to_be_retained[chromosome]=[]
        intron_motif,uniq_reads,mm_reads=int(intron_motif),int(uniq_reads),int(mm_reads)
        
        if prev_chromosome=="":
            prev_chromosome=chromosome
            prev_junction=j_start+"-"+j_end
        if j_start+"-"+j_end == prev_junction and prev_chromosome==chromosome:
            same_junction.append([intron_motif,uniq_reads,mm_reads,strand])
        else:
            # Non-cannonical junctions are retained only if they have very high read support in at least 3 samples of the SAME project
            try :        
                min_uniq_reads=7 if same_junction[0][0]==0 else 3
                min_mm_reads=10 if same_junction[0][0]==0 else 6
                if len([row for row in same_junction if row[1]<min_uniq_reads and row[2]<min_mm_reads])>0:
                    if num_samples<=1:
                        junctions_to_be_discarded[chromosome].append(prev_junction)                        
                    else:
                        """if condition=="dummylight" and prev_junction=="10042104-10044947" and chromosome=="1":
                            print(num_samples<=3,len([row for row in same_junction if row[1]>=min_uniq_reads or row[2]>=min_mm_reads]),[row for row in same_junction if row[1]>=min_uniq_reads or row[2]>=min_mm_reads])"""
                        if num_samples<=3 and len([row for row in same_junction if row[1]>=min_uniq_reads or row[2]>=min_mm_reads])==0:
                            """if condition=="dummylight" and prev_junction=="10042104-10044947" and chromosome=="1":print("INSIDE")"""
                            junctions_to_be_discarded[chromosome].append(prev_junction)
                        elif num_samples>3:
                            if len(same_junction)<3:
                                junctions_to_be_discarded[chromosome].append(prev_junction)
                        else:
                            #junctions_to_be_retained[chromosome].append(prev_junction+"_"+strand)
                            junctions_to_be_retained[chromosome].append(prev_line)
                else:
                    #junctions_to_be_retained[chromosome].append(prev_junction+"_"+strand)
                    junctions_to_be_retained[chromosome].append(prev_line)
            except IndexError:
                print("ERROR",condition,chromosome,prev_junction,same_junction)
                return [],[]
            same_junction=[]
            same_junction.append([intron_motif,uniq_reads,mm_reads])
            if prev_chromosome!=chromosome:
                prev_chromosome=chromosome
        prev_junction=j_start+"-"+j_end
        prev_line=line
    fhr.close()
    return condition,junctions_to_be_discarded,junctions_to_be_retained
        
def writeNewSJDBFileInParallel(eachinput):
    sjdbout_filename_merged_condition,new_sjdbout_filename_merged_condition,junctions_to_be_retained=eachinput
    fhr=open(sjdbout_filename_merged_condition,"r")
    fhw=open(new_sjdbout_filename_merged_condition,"w")
    for chromosome in junctions_to_be_retained:
        for line in junctions_to_be_retained[chromosome]:
            fhw.write(line)
    fhw.close()
    
    cmd="mv "+new_sjdbout_filename_merged_condition+" "+sjdbout_filename_merged_condition
    os.system(cmd)

def selectHighConfidenceSpliceJunctions(options,round,condition):
    condition,junctions_to_be_discarded,junctions_to_be_retained=selectHighConfidenceSpliceJunctionsPerCondition([condition,options,round])
    sjdbout_filename_merged_condition=options.output_star+"/"+condition+"_round"+str(round)+"_SJ.out.tab"
    new_sjdbout_filename_merged_condition=options.output_star+"/"+condition+"_round"+str(round)+"_SJ.out.tab.new"
    writeNewSJDBFileInParallel([sjdbout_filename_merged_condition,new_sjdbout_filename_merged_condition,junctions_to_be_retained])

def combineSpliceJunctionDatabases(options,rounds,condition):
    
    cmd="cat "
    for round in rounds:
        cmd+=options.output_star+"/"+condition+"_round"+str(round)+"_SJ.out.tab "
    if len(rounds)==2:
        cmd+=" > "+options.output_star+"/"+condition+"_round1_and_round2_SJ.out.tab "
    elif len(rounds)==3:
        cmd+=" > "+options.output_star+"/"+condition+"_round1_and_round2_and_round3_SJ.out.tab "
    elif len(rounds)==4:
        cmd+=" > "+options.output_star+"/"+condition+"_round1_and_round2_and_round3_and_round4_SJ.out.tab "
    os.system(cmd)

def checkMappingRateToDiscardPoorRuns(options,condition,Run,logging_mutex,logger_proxy):
    """
    """
    total_reads1,umr1,mmr1=pullOutMappingInformation(options.output_star+"/"+Run+"_round1_Log.final.out")
    total_reads2,umr2,mmr2=pullOutMappingInformation(options.output_star+"/"+Run+"_round2_Log.final.out")
    if total_reads2>0:
        total_mapped_reads_percentage=(umr2+mmr2)/total_reads2
    else:
        total_mapped_reads_percentage=0
    with logging_mutex:
        logger_proxy.info("Mapping rate in round2 "+condition+" "+Run+" "+str(total_mapped_reads_percentage))
    if total_mapped_reads_percentage<0.2:
        return 0
    else:
        return 1

def mergeAllAlignments(options,condition):
    pool = multiprocessing.Pool(processes=int(options.cpu))
    allinputs=[]
    
    for Run in options.mrna_md[condition]:
        bam_filename_round1=options.output_star+"/"+Run+"_round1_Aligned.sortedByCoord.out.bam"
        bam_filename_round2=options.output_star+"/"+Run+"_round2_Aligned.sortedByCoord.out.bam"
        bam_filename_round3=options.output_star+"/"+Run+"_round3_Aligned.sortedByCoord.out.bam"
        #bam_filename_round4=options.output_star+"/"+Run+"_round4_Aligned_portcullis_filtered.sortedByCoord.out.bam"
        bam_filename_round4=options.output_star+"/"+Run+"_round4_Aligned.sortedByCoord.out.bam"
        bam_filename_round5=options.output_star+"/"+Run+"_olego_round5.sorted.bam"
        bam_filename_final=options.output_star+"/"+Run+"_final.sortedByCoord.out.bam"
        
        if os.path.exists(bam_filename_final)==True:
            continue
        cmd="samtools "+" merge -@ "+str(options.cpu)+" "
        cmd+=" -f "+bam_filename_final+" "
        cmd+=bam_filename_round1+" "
        cmd+=bam_filename_round2+" "
        cmd+=bam_filename_round3+" "
        cmd+=bam_filename_round4+" "
        cmd+=bam_filename_round5+" "
        allinputs.append([Run,cmd])
    if len(allinputs)>0:
        pool.map(runCommand,allinputs)
            
    for Run in options.mrna_md[condition]:
        bam_filename_final=options.output_star+"/"+Run+"_final.sortedByCoord.out.bam"  
        cmd="samtools "+" index -@ "+str(options.cpu)+" "
        cmd+=bam_filename_final
        if os.path.exists(bam_filename_final+".bai")==False:
            os.system(cmd)
        cmd="samtools "+" index -c -@ "+str(options.cpu)+" "
        cmd+=bam_filename_final
        if os.path.exists(bam_filename_final+".csi")==False:
            os.system(cmd)

def performRandomSelectionOfReads(proportion_of_reads_to_be_selected,input_filename1,output_filename1,input_filename2=None,output_filename2=None):
    if input_filename2==None and output_filename2==None: # single ended
        fhr1 = open(input_filename1,"r")
        fhw1 = open(output_filename1,"w")
        while True:
            line1_1 = fhr1.readline()
            if not line1_1:break
            line1_2 = fhr1.readline()
            line1_3 = fhr1.readline()
            line1_4 = fhr1.readline()
            if random.uniform(0,1) < proportion_of_reads_to_be_selected:
                fhw1.write(line1_1)
                fhw1.write(line1_2)
                fhw1.write(line1_3)
                fhw1.write(line1_4)
        fhr1.close()
        fhw1.close()
    else:
        fhr1 = open(input_filename1,"r")
        fhr2 = open(input_filename2,"r")
        fhw1 = open(output_filename1,"w")
        fhw2 = open(output_filename2,"w")
        while True:
            line1_1 = fhr1.readline()
            if not line1_1:break
            line1_2 = fhr1.readline()
            line1_3 = fhr1.readline()
            line1_4 = fhr1.readline()
            
            line2_1 = fhr2.readline()
            line2_2 = fhr2.readline()
            line2_3 = fhr2.readline()
            line2_4 = fhr2.readline()
            if random.uniform(0,1) < proportion_of_reads_to_be_selected:
                fhw1.write(line1_1)
                fhw1.write(line1_2)
                fhw1.write(line1_3)
                fhw1.write(line1_4)
                
                fhw2.write(line2_1)
                fhw2.write(line2_2)
                fhw2.write(line2_3)
                fhw2.write(line2_4)
        fhr1.close()
        fhw1.close()
        fhr2.close()
        fhw2.close()

def selectReadsAtRandom(sraids_filenames,location_directory,options,condition,logger_proxy,logging_mutex):
    """
    Randomly select reads from fastq
    """
    proportion_of_reads_to_be_selected=0.1
    fhr = open(sraids_filenames,"r")
    for Run in fhr:
        Run=Run.strip()
        with logging_mutex:
            logger_proxy.info(f"Randomly selecting {proportion_of_reads_to_be_selected} reads from {Run}")
        if options.mrna_md[condition][Run]["Ended"]=="SE":
            input_filename1 = location_directory+"/"+Run+".fastq"
            output_filename1 = location_directory+"/"+Run+".fastq.reduced"
            performRandomSelectionOfReads(proportion_of_reads_to_be_selected,input_filename1,output_filename1)
            os.system(f"mv {output_filename1} {input_filename1}")
        elif options.mrna_md[condition][Run]["Ended"]=="PE":
            input_filename1 = location_directory+"/"+Run+"_1.fastq"
            output_filename1 = location_directory+"/"+Run+"_1.fastq.reduced"
            input_filename2 = location_directory+"/"+Run+"_2.fastq"
            output_filename2 = location_directory+"/"+Run+"_2.fastq.reduced"
            performRandomSelectionOfReads(proportion_of_reads_to_be_selected,input_filename1,output_filename1,input_filename2,output_filename2)
            os.system(f"mv {output_filename1} {input_filename1}")
            os.system(f"mv {output_filename2} {input_filename2}")
    fhr.close()

def alignReadsAndMergeOutput(options,logger_proxy,logging_mutex):
    for condition_num,condition in enumerate(options.mrna_md):
        with logging_mutex:
            logger_proxy.info("Started processing data for "+condition)
        
        flag=0
        for Run in options.mrna_md[condition]:
            if os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==False and os.path.exists(options.output_star+"/"+Run+"_for_psiclass.bam")==False:
                flag=1
        if flag==0:continue
        #########################################################################################################
        # Download missing data from NCBI
        #########################################################################################################
        fhw=open(options.temp_dir+"/download_these_runs","w")
        download_these=[]
        for Run in options.mrna_md[condition]:
            if os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True and os.path.exists(options.output_star+"/"+Run+"_round4_Aligned.sortedByCoord.out.bam")==True:continue
            if options.mrna_md[condition][Run]["downloaded_from_NCBI"]==0:continue
            #print(options.mrna_md[condition][Run])
            if options.mrna_md[condition][Run]["Ended"]=="SE":
                if os.path.exists(options.raw_data_downloaded_from_NCBI+"/"+Run+".fastq")==False and os.path.exists(options.output_fasta_N_removed+"/"+Run+".fasta")==False and os.path.exists(options.raw_data_downloaded_from_NCBI+"/"+Run+".fq")==False:
                    fhw.write(Run+"\n")
                    download_these.append(Run)
                    #options.mrna_md[condition][Run]["location_directory"]=options.raw_data_downloaded_from_NCBI
            elif options.mrna_md[condition][Run]["Ended"]=="PE":
                if (os.path.exists(options.raw_data_downloaded_from_NCBI+"/"+Run+"_1.fastq")==False or os.path.exists(options.raw_data_downloaded_from_NCBI+"/"+Run+"_2.fastq")==False) \
                and (os.path.exists(options.output_fasta_N_removed+"/"+Run+"_1.fasta")==False or os.path.exists(options.output_fasta_N_removed+"/"+Run+"_2.fasta")==False) \
                and (os.path.exists(options.raw_data_downloaded_from_NCBI+"/"+Run+"_1.fq")==False or os.path.exists(options.raw_data_downloaded_from_NCBI+"/"+Run+"_2.fq")==False):
                    fhw.write(Run+"\n")
                    download_these.append(Run)
                    #options.mrna_md[condition][Run]["location_directory"]=options.raw_data_downloaded_from_NCBI
        fhw.close()
        #print(open(options.temp_dir+"/download_these_runs","r").read())
        with logging_mutex:
            logger_proxy.info("Downloading missing data from NCBI started")
        
        if len(download_these)>0:
            cmd=options.softwares["download_and_dump_fastq_from_SRA"]
            cmd+=" -s "+options.temp_dir+"/download_these_runs "
            cmd+=" -o "+options.raw_data_downloaded_from_NCBI
            cmd+=" -n "+str(options.cpu)
            cmd+=" > "+options.temp_dir+"/download_these_runs.output "
            cmd+=" 2> "+options.temp_dir+"/download_these_runs.error "
            
            with logging_mutex:
                logger_proxy.info("Running command - "+cmd)
            os.system(cmd)
        with logging_mutex:
            logger_proxy.info("Downloading missing data from NCBI finished for "+condition)
        
        
        if options.run_tests==True:
            selectReadsAtRandom(options.temp_dir+"/download_these_runs",options.raw_data_downloaded_from_NCBI,options,condition,logger_proxy,logging_mutex)
        rerun_selection_of_high_confidence_SJ=0 # Assumes that all the Runs have been aligned
        #########################################################################################################
        # Align reads with STAR round1 & compress adapter trimmed reads
        #########################################################################################################
        for Run in options.mrna_md[condition]:
            if (os.path.exists(options.output_star+"/"+Run+"_round1_Aligned.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_round1_Aligned.sortedByCoord.out.bam",options)==0) or os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:continue
            rerun_selection_of_high_confidence_SJ=1
            alignReadsWithSTARRound1(options,Run,options.mrna_md[condition][Run]["Ended"],condition,logger_proxy,logging_mutex)
            
        with logging_mutex:
            logger_proxy.info("Mapping of reads for round1 completed for "+condition)
        
        if rerun_selection_of_high_confidence_SJ==1:
            selectHighConfidenceSpliceJunctions(options,1,condition)
            with logging_mutex:
                logger_proxy.info("Selecting high confidence junctions after round1 mapping completed for "+condition)
        else:
            with logging_mutex:
                logger_proxy.info("High confidence junctions after round1 mapping is present for "+condition)
         
        #########################################################################################################
        # Remove files to clean up space
        #########################################################################################################
        if options.no_cleanup==False:
            # Remove all raw fastq files downloaded from NCBI
            #options.space_saved+=sum(f.stat().st_size for f in Path(options.raw_data_downloaded_from_NCBI).glob('**/*') if f.is_file() )
            cmd="rm -f "+options.raw_data_downloaded_from_NCBI+"/*"
            os.system(cmd)
            
            # Remove all Rcorrected read files
            #options.space_saved+=sum(f.stat().st_size for f in Path(options.error_corrected_raw_data).glob('**/*') if f.is_file() )
            cmd="rm -f "+options.error_corrected_raw_data+"/*"
            os.system(cmd)
            with logging_mutex:
                logger_proxy.info("Raw read download from NCBI cleanup completed for "+condition)
        
        #########################################################################################################
        # Align reads with STAR round2
        #########################################################################################################
        for Run in options.mrna_md[condition]:
            if (os.path.exists(options.output_star+"/"+Run+"_round2_Aligned.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_round2_Aligned.sortedByCoord.out.bam",options)==0) or os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:continue
            rerun_selection_of_high_confidence_SJ=1
            alignReadsWithSTARRound2(options,Run,options.mrna_md[condition][Run]["Ended"],condition,logger_proxy,logging_mutex)
        with logging_mutex:
            logger_proxy.info("Mapping of reads for round2 completed for "+condition)
        
        if rerun_selection_of_high_confidence_SJ==1:
            selectHighConfidenceSpliceJunctions(options,2,condition)
            combineSpliceJunctionDatabases(options,[1,2],condition)
            with logging_mutex:
                logger_proxy.info("Selecting high confidence junctions after round2 mapping completed for "+condition)
        else:
            with logging_mutex:
                logger_proxy.info("High confidence junctions after round2 mapping is present for "+condition)
                
        #########################################################################################################
        # Flag runs with very low number of reads mapping
        #########################################################################################################
        default_align_these_runs=[]
        for Run in options.mrna_md[condition]:
            if os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:continue
            if checkMappingRateToDiscardPoorRuns(options,condition,Run,logging_mutex,logger_proxy)==0:
                default_align_these_runs.append(Run)
        
        with logging_mutex:
            logger_proxy.info("Resorting to alignment with relaxed parameters for these runs due to poor mapping "+",".join(default_align_these_runs))
        
        for Run in default_align_these_runs:
            if os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:continue
            alignReadsWithSTARRelaxed(options, Run, options.mrna_md[condition][Run]["Ended"], condition, logger_proxy, logging_mutex)
        
        #########################################################################################################
        # Align reads with STAR round3
        #########################################################################################################
        for Run in options.mrna_md[condition]:
            if (os.path.exists(options.output_star+"/"+Run+"_round3_Aligned.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_round3_Aligned.sortedByCoord.out.bam",options)==0) or os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:continue
            if Run in default_align_these_runs:continue
            alignReadsWithSTARRound3(options,Run,options.mrna_md[condition][Run]["Ended"],condition,logger_proxy,logging_mutex)
            rerun_selection_of_high_confidence_SJ=1
        with logging_mutex:
            logger_proxy.info("Mapping of reads for round3 completed for "+condition)
            
        if rerun_selection_of_high_confidence_SJ==1:
            selectHighConfidenceSpliceJunctions(options,3,condition)
            combineSpliceJunctionDatabases(options,[1,2,3],condition)
            with logging_mutex:
                logger_proxy.info("Selecting high confidence junctions after round3 mapping completed for "+condition)
        else:
            with logging_mutex:
                logger_proxy.info("High confidence junctions after round3 mapping is present for "+condition)
        
        #########################################################################################################
        # Align reads with STAR round4
        #########################################################################################################
        for Run in options.mrna_md[condition]:
            if Run in default_align_these_runs:continue
            if (os.path.exists(options.output_star+"/"+Run+"_round4_Aligned.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_round4_Aligned.sortedByCoord.out.bam",options)==0):continue
            alignReadsWithSTARRound4(options,Run,options.mrna_md[condition][Run]["Ended"],condition,logger_proxy,logging_mutex)
            rerun_selection_of_high_confidence_SJ=1
        with logging_mutex:
            logger_proxy.info("Mapping of reads for round4 completed for "+condition)
            
        if rerun_selection_of_high_confidence_SJ==1:
            selectHighConfidenceSpliceJunctions(options,4,condition)
            combineSpliceJunctionDatabases(options,[1,2,3,4],condition)
            with logging_mutex:
                logger_proxy.info("Selecting high confidence junctions after round4 mapping completed for "+condition)
        else:
            with logging_mutex:
                logger_proxy.info("High confidence junctions after round4 mapping is present for "+condition)
        
        #########################################################################################################
        # Align reads with OLego round5
        #########################################################################################################
        flag=1
        for Run_num,Run in enumerate(options.mrna_md[condition]):
            if Run in default_align_these_runs:continue
            fhr=open(options.output_star+"/"+Run+"_round4_Log.final.out","r")
            for line in fhr:
                if "Number of reads unmapped: too short" in line:
                    num_unmapped_reads=int(line.strip().split()[-1])
            fhr.close()
            if num_unmapped_reads>5000000:
                # Create dummy olego file to assist in merge operation
                os.system("samtools view -Hb "+options.output_star+"/"+Run+"_round4_Aligned.sortedByCoord.out.bam "+" >  "+options.output_star+"/"+Run+"_olego_round5.sorted.bam")
                with logging_mutex:
                    logger_proxy.info("Mapping with OLego for micro-exon detection skipped for "+Run+" due to huge number of reads")
            else:
                if flag==1:
                    start_time=time.time()
                    flag=alignReadsWithOLegoRound5(options,Run,options.mrna_md[condition][Run]["Ended"],condition,logger_proxy,logging_mutex,2)
                    end_time=time.time()
                else:
                    if end_time-start_time<7200:
                        alignReadsWithOLegoRound5(options,Run,options.mrna_md[condition][Run]["Ended"],condition,logger_proxy,logging_mutex,2)
                    else:
                        alignReadsWithOLegoRound5(options,Run,options.mrna_md[condition][Run]["Ended"],condition,logger_proxy,logging_mutex,4)
        with logging_mutex:
            logger_proxy.info("Mapping with OLego for micro-exon detection completed for "+condition)
        
        #########################################################################################################
        # Merge alignments from each of the 5 rounds
        #########################################################################################################
        mergeAllAlignments(options,condition)
        with logging_mutex:
            logger_proxy.info("Merging of alignments from all rounds of mapping completed for "+condition)
        
        #########################################################################################################
        # Remove all alignments from all the rounds to save space
        #########################################################################################################
        if options.no_cleanup==False:
            for Run in options.mrna_md[condition]:
                cmd="rm -rf "
                cmd+=options.output_star+"/"+Run+"_round1_Aligned.sortedByCoord.out.bam"
                os.system(cmd)
                
                cmd="rm -rf "
                cmd+=options.output_star+"/"+Run+"_round2_Aligned.sortedByCoord.out.bam"
                os.system(cmd)
                
                cmd="rm -rf "
                cmd+=options.output_star+"/"+Run+"_round3_Aligned.sortedByCoord.out.bam"
                os.system(cmd)
                                
                cmd="rm -rf "
                cmd+=options.output_star+"/"+Run+"_round*__STARgenome"
                os.system(cmd)
                
                cmd="rm -rf "
                cmd+=options.output_star+"/"+Run+"_round*_Unmapped.out.mate*"
                os.system(cmd)
                
                cmd="rm -rf "
                cmd+=options.output_star+"/"+Run+"_olego_round5.sorted.bam"
                os.system(cmd)
                    
                cmd="rm -rf "
                cmd+=options.output_star+"/"+Run+"*STARtmp"
                os.system(cmd)
                
            with logging_mutex:
                logger_proxy.info("Removing intermediate alignment files completed for "+condition)
                
        if options.preserve_raw_input_data==False:
            for Run in options.mrna_md[condition]:
                cmd="rm -rf "+options.output_fasta_N_removed+"/"+Run+"*"
                os.system(cmd)
        
        #########################################################################################################
        # Log final completion for the condition
        #########################################################################################################
        with logging_mutex:
            logger_proxy.info("Mapping of all runs completed for "+condition)
