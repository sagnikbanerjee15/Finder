

from scripts.fileReadWriteOperations import *
from scripts.runCommand import *
import os
import sys

def alignReadsWithSTARRound1(options,Run,ended,condition,logger_proxy,logging_mutex):
    if (os.path.exists(options.output_star+"/"+Run+"_round1_Aligned.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_round1_Aligned.sortedByCoord.out.bam",options)==0) or os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:return
    cmd="STAR "
    cmd+=" --runThreadN "+str(options.cpu)
    cmd+=" --genomeDir "+options.genome_dir_star
    cmd+=" --outSAMtype BAM SortedByCoordinate "
    cmd+=" --outFilterMultimapNmax 500 " 
    cmd+=" --outFilterMismatchNmax 0  " # For round 1 no mismatches are allowed --> done to capture the most confident junction
    cmd+=" --alignIntronMin 20  "
    cmd+=" --alignIntronMax 10000 "
    cmd+=" --limitBAMsortRAM 107374182400"
    cmd+=" --alignEndsType EndToEnd " # Disallow soft clipping 
    cmd+=" --outSAMprimaryFlag AllBestScore "
    #cmd+=" --outFilterScoreMinOverLread 0.90 "
    #cmd+=" --outFilterMatchNminOverLread 0.90 "
    cmd+=" --outSJfilterOverhangMin 12 12 12 12 "
    cmd+=" --outSAMattributes NH HI AS nM NM MD jM jI XS "
    cmd+=" --outReadsUnmapped Fastx "
    if options.star_shared_mem == True:
        cmd+=" --genomeLoad LoadAndKeep "
    cmd+=" --outFileNamePrefix "+options.output_star+"/"+Run+"_round1_"
    cmd+=" --outSJfilterCountUniqueMin 1 1 1 1 "
    cmd+=" --outSJfilterCountTotalMin 1 1 1 1 "
    cmd+=" --alignSJoverhangMin 12 " # For round 1 a high overhang ensures reads mapping with high confidence    
    cmd+=" --outSAMattrRGline ID:1 "
    
    if ended=="SE":
        cmd+=" --readFilesIn "+options.mrna_md[condition][Run]["location_directory"]+"/"+Run+".fastq"
    else:
        cmd+=" --readFilesIn "+options.mrna_md[condition][Run]["location_directory"]+"/"+Run+"_1.fastq "
        cmd+=" "+options.mrna_md[condition][Run]["location_directory"]+"/"+Run+"_2.fastq "
    cmd+=" > "+options.output_star+"/"+Run+"_round1.output"
    cmd+=" 2> "+options.output_star+"/"+Run+"_round1.error"
    os.system(cmd)
    
    cmd="rm "
    cmd+=options.output_star+"/"+Run+"_round1_Log.progress.out "
    cmd+=options.output_star+"/"+Run+"_round1_Log.out "
    os.system(cmd)
    
    if options.star_shared_mem == True:
        # Remove a pre-loaded genome
        cmd  = "STAR "
        cmd += f" --runThreadN {options.cpu} "
        cmd += f" --genomeLoad Remove "
        cmd += f" --genomeDir {options.genome_dir_star} "
        os.system(cmd)
    
    with logging_mutex:
        logger_proxy.info("STAR Round1 run for "+Run+" completed")

def alignReadsWithSTARRelaxed(options,Run,ended,condition,logger_proxy,logging_mutex):
    
    if os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam",options)==0:return
    cmd="STAR "
    cmd+=" --runThreadN "+str(options.cpu)
    cmd+=" --genomeDir "+options.genome_dir_star
    """if Run_num==0:
        cmd+=" --genomeDir "+options.genome_dir_star
    else:
        cmd+=" --genomeDir "+options.output_star+"/"+bioproject_to_run[bioproject][0]+"_round2__STARgenome"
    """
    cmd+=" --outSAMtype BAM SortedByCoordinate "
    cmd+=" --outFilterMultimapNmax 500 " 
    cmd+=" --outFilterMismatchNmax 5  " # For relaxed 5 mismatches are allowed 
    cmd+=" --alignIntronMin 20  "
    cmd+=" --alignIntronMax 10000 "
    
    cmd+=" --limitBAMsortRAM 107374182400 "
    cmd+=" --alignEndsType Local " # ALLOWS soft clipping 
    cmd+=" --outSAMprimaryFlag AllBestScore "
    cmd+=" --outFilterScoreMinOverLread 0.60"
    cmd+=" --outFilterMatchNminOverLread 0.60 "
    cmd+=" --outSJfilterOverhangMin 8 8 8 8 "
    cmd+=" --outSAMattributes NH HI AS nM NM MD jM jI XS "
    cmd+=" --outReadsUnmapped Fastx "
    cmd+=" --outSAMattrRGline ID:6 "
    """if Run_num==0:
        cmd+=" --sjdbInsertSave All "
        cmd+=" --sjdbFileChrStartEnd "+options.output_star+"/"+condition+"_round1_SJ.out.tab "
    """
    
    """if Run_num!=0:
        if options.star_shared_mem == True:
            cmd+=" --genomeLoad LoadAndKeep "
    """    
    cmd+=" --outFileNamePrefix "+options.output_star+"/"+Run+"_final_"
    cmd+=" --outSJfilterCountUniqueMin 1 1 1 1 " 
    cmd+=" --outSJfilterCountTotalMin 1 1 1 1 " 
    cmd+=" --alignSJoverhangMin 8 "
    
    if ended=="SE":
        cmd+=" --readFilesIn "+options.mrna_md[condition][Run]["location_directory"]+"/"+Run+".fastq"
    else:
        cmd+=" --readFilesIn "+options.mrna_md[condition][Run]["location_directory"]+"/"+Run+"_1.fastq "
        cmd+=" "+options.mrna_md[condition][Run]["location_directory"]+"/"+Run+"_2.fastq "
    cmd+=" > "+options.output_star+"/"+Run+"_relaxed.output"
    cmd+=" 2> "+options.output_star+"/"+Run+"_relaxed.error"
    #print(cmd)
    #sys.stdout.flush()
    os.system(cmd)
    
    cmd="rm "
    cmd+=options.output_star+"/"+Run+"_final_Log.progress.out "+options.output_star+"/"+Run+"_final_Log.out"
    os.system(cmd)
    
    cmd="mv "
    cmd+=options.output_star+"/"+Run+"_final_Aligned.sortedByCoord.out.bam "
    cmd+=options.output_star+"/"+Run+"_final.sortedByCoord.out.bam"
    os.system(cmd)
    
    cmd="mv "
    cmd+=options.output_star+"/"+Run+"_final_Unmapped.out.mate1 "
    cmd+=options.output_star+"/"+Run+"_round3_Unmapped.out.mate1"
    os.system(cmd)
    if ended=="PE":
        cmd="mv "
        cmd+=options.output_star+"/"+Run+"_final_Unmapped.out.mate2 "
        cmd+=options.output_star+"/"+Run+"_round3_Unmapped.out.mate2"
        os.system(cmd)
    
    cmd="mv "
    cmd+=options.output_star+"/"+Run+"_final_Log.final.out "
    cmd+=options.output_star+"/"+Run+"_round3_Log.final.out"
    os.system(cmd)
    
    if options.star_shared_mem == True:
        # Remove a pre-loaded genome
        cmd  = "STAR "
        cmd += f" --runThreadN {options.cpu} "
        cmd += f" --genomeLoad Remove "
        cmd += f" --genomeDir {options.genome_dir_star} "
        os.system(cmd)

    with logging_mutex:
        logger_proxy.info("STAR relaxed alignment run for "+Run+" completed")
              
def alignReadsWithSTARRound2(options,Run,ended,condition,logger_proxy,logging_mutex):
    
    if (os.path.exists(options.output_star+"/"+Run+"_round2_Aligned.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_round2_Aligned.sortedByCoord.out.bam",options)==0) or os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:return
    cmd="STAR "
    cmd+=" --runThreadN "+str(options.cpu)
    cmd+=" --genomeDir "+options.genome_dir_star
    """if Run_num==0:
        cmd+=" --genomeDir "+options.genome_dir_star
    else:
        cmd+=" --genomeDir "+options.output_star+"/"+bioproject_to_run[bioproject][0]+"_round2__STARgenome"
    """
    cmd+=" --outSAMtype BAM SortedByCoordinate "
    cmd+=" --outFilterMultimapNmax 500 " 
    cmd+=" --outFilterMismatchNmax 2  " # For round 2 TWO mismatches are allowed 
    cmd+=" --alignIntronMin 20  "
    cmd+=" --alignIntronMax 10000 "
    cmd+=" --limitBAMsortRAM 107374182400 "
    cmd+=" --alignEndsType Local " # ALLOWS soft clipping 
    cmd+=" --outSAMprimaryFlag AllBestScore "
    cmd+=" --outFilterScoreMinOverLread 0.95 "
    cmd+=" --outFilterMatchNminOverLread 0.95 "
    cmd+=" --outSJfilterOverhangMin 8 8 8 8 "
    cmd+=" --outSAMattributes NH HI AS nM NM MD jM jI XS "
    cmd+=" --outReadsUnmapped Fastx "
    cmd+=" --outSAMattrRGline ID:2 "
    """if Run_num==0:
        cmd+=" --sjdbInsertSave All "
        cmd+=" --sjdbFileChrStartEnd "+options.output_star+"/"+condition+"_round1_SJ.out.tab "
    """
    cmd+=" --sjdbFileChrStartEnd "+options.output_star+"/"+condition+"_round1_SJ.out.tab "
    
    """if Run_num!=0:
        if options.star_shared_mem == True:
            cmd+=" --genomeLoad LoadAndKeep "
    """    
    cmd+=" --outFileNamePrefix "+options.output_star+"/"+Run+"_round2_"
    cmd+=" --outSJfilterCountUniqueMin 1 1 1 1 " 
    cmd+=" --outSJfilterCountTotalMin 1 1 1 1 " 
    cmd+=" --alignSJoverhangMin 12 " # For round 2 same high overhang is used as in round 1
    cmd+=" --alignSJDBoverhangMin 8 " # 
    
    if ended=="SE":
        cmd+=" --readFilesIn "+options.output_star+"/"+Run+"_round1_Unmapped.out.mate1"
    else:
        cmd+=" --readFilesIn "+options.output_star+"/"+Run+"_round1_Unmapped.out.mate1 "
        cmd+=" "+options.output_star+"/"+Run+"_round1_Unmapped.out.mate2"
    cmd+=" > "+options.output_star+"/"+Run+"_round2.output"
    cmd+=" 2> "+options.output_star+"/"+Run+"_round2.error"
    #print(cmd)
    #sys.stdout.flush()
    os.system(cmd)
    
    if options.star_shared_mem == True:
        # Remove a pre-loaded genome
        cmd  = "STAR "
        cmd += f" --runThreadN {options.cpu} "
        cmd += f" --genomeLoad Remove "
        cmd += f" --genomeDir {options.genome_dir_star} "
        os.system(cmd)
    
    cmd="rm "
    cmd+=options.output_star+"/"+Run+"_round2_Log.progress.out "+options.output_star+"/"+Run+"_round2_Log.out"
    os.system(cmd)
    
    with logging_mutex:
        logger_proxy.info("STAR Round2 run for "+Run+" completed")
    
def alignReadsWithSTARRound3(options,Run,ended,condition,logger_proxy,logging_mutex):
    
    if (os.path.exists(options.output_star+"/"+Run+"_round3_Aligned.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_round3_Aligned.sortedByCoord.out.bam",options)==0) or os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:return
    cmd="STAR "
    cmd+=" --runThreadN "+str(options.cpu)
    cmd+=" --genomeDir "+options.genome_dir_star
    """if Run_num==0:
        cmd+=" --genomeDir "+options.genome_dir_star
    else:
        cmd+=" --genomeDir "+options.output_star+"/"+bioproject_to_run[bioproject][0]+"_round2__STARgenome"
    """
    cmd+=" --outSAMtype BAM SortedByCoordinate "
    cmd+=" --outFilterMultimapNmax 500 " 
    cmd+=" --outFilterMismatchNmax 3  " # For round 3 THREE mismatches are allowed 
    cmd+=" --alignIntronMin 20  "
    cmd+=" --alignIntronMax 10000 "
    cmd+=" --limitBAMsortRAM 107374182400"
    cmd+=" --alignEndsType Local " # ALLOWS soft clipping 
    cmd+=" --outSAMprimaryFlag AllBestScore "
    cmd+=" --outFilterScoreMinOverLread 0.95 "
    cmd+=" --outFilterMatchNminOverLread 0.95 "
    cmd+=" --outSJfilterOverhangMin 8 8 8 8 "
    cmd+=" --outSAMattributes NH HI AS nM NM MD jM jI XS "
    cmd+=" --outReadsUnmapped Fastx "
    cmd+=" --outSAMattrRGline ID:3 "
    """if Run_num==0:
        cmd+=" --sjdbInsertSave All "
        cmd+=" --sjdbFileChrStartEnd "+options.output_star+"/"+condition+"_round2_SJ.out.tab "
    """
    cmd+=" --sjdbFileChrStartEnd "+options.output_star+"/"+condition+"_round2_SJ.out.tab "
    
    """if Run_num!=0:
        if options.star_shared_mem == True:
            cmd+=" --genomeLoad LoadAndKeep "
    """    
    cmd+=" --outFileNamePrefix "+options.output_star+"/"+Run+"_round3_"
    cmd+=" --outSJfilterCountUniqueMin 1 1 1 1 " 
    cmd+=" --outSJfilterCountTotalMin 1 1 1 1 " 
    cmd+=" --alignSJoverhangMin 1000 " # Very high theshold used to prevent creation of new junctions
    cmd+=" --alignSJDBoverhangMin 8 " # 
    
    if ended=="SE":
        cmd+=" --readFilesIn "+options.output_star+"/"+Run+"_round2_Unmapped.out.mate1"
    else:
        cmd+=" --readFilesIn "+options.output_star+"/"+Run+"_round2_Unmapped.out.mate1 "
        cmd+=" "+options.output_star+"/"+Run+"_round2_Unmapped.out.mate2"
    cmd+=" > "+options.output_star+"/"+Run+"_round3.output"
    cmd+=" 2> "+options.output_star+"/"+Run+"_round3.error"
    #print(cmd)
    sys.stdout.flush()
    os.system(cmd)
    
    if options.star_shared_mem == True:
        # Remove a pre-loaded genome
        cmd  = "STAR "
        cmd += f" --runThreadN {options.cpu} "
        cmd += f" --genomeLoad Remove "
        cmd += f" --genomeDir {options.genome_dir_star} "
        os.system(cmd)
    
    cmd="rm "
    cmd+=options.output_star+"/"+Run+"_round3_Log.progress.out "+options.output_star+"/"+Run+"_round3_Log.out"
    os.system(cmd)
    
    with logging_mutex:
        logger_proxy.info("STAR Round3 run for "+Run+" completed")

def alignReadsWithSTARRound4(options,Run,ended,condition,logger_proxy,logging_mutex): 
    if (os.path.exists(options.output_star+"/"+Run+"_round4_Aligned.sortedByCoord.out.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_round4_Aligned.sortedByCoord.out.bam",options)==0):return
    cmd="STAR "
    cmd+=" --runThreadN "+str(options.cpu)
    cmd+=" --genomeDir "+options.genome_dir_star
    cmd+=" --outSAMtype BAM SortedByCoordinate "
    cmd+=" --outFilterMultimapNmax 500 " 
    cmd+=" --outFilterMismatchNmax 2  " # For round 4 TWO mismatches are allowed 
    cmd+=" --alignIntronMin 10000  "
    cmd+=" --alignIntronMax 100000000 "
    cmd+=" --limitBAMsortRAM 107374182400"
    cmd+=" --alignEndsType Local " # ALLOWS soft clipping 
    cmd+=" --outSAMprimaryFlag AllBestScore "
    cmd+=" --outFilterScoreMinOverLread 0.95 "
    cmd+=" --outFilterMatchNminOverLread 0.95 "
    cmd+=" --outSJfilterOverhangMin 8 8 8 8 "
    cmd+=" --outSAMattributes NH HI AS nM NM MD jM jI XS "
    cmd+=" --outReadsUnmapped Fastx "
    cmd+=" --outSAMattrRGline ID:4 "
    if options.star_shared_mem == True:
        cmd+=" --genomeLoad LoadAndKeep "
    cmd+=" --outFileNamePrefix "+options.output_star+"/"+Run+"_round4_"
    cmd+=" --outSJfilterCountUniqueMin 3 3 3 3 " # Stricter bounds than round1
    cmd+=" --outSJfilterCountTotalMin 5 5 5 5 " # Stricter bounds than round1
    cmd+=" --alignSJoverhangMin 8 " # For round 4 same high overhang is used as in round 1
    #cmd+=" --alignSJDBoverhangMin 8 "
    
    if ended=="SE":
        cmd+=" --readFilesIn "+options.output_star+"/"+Run+"_round3_Unmapped.out.mate1"
    else:
        cmd+=" --readFilesIn "+options.output_star+"/"+Run+"_round3_Unmapped.out.mate1"
        cmd+=" "+options.output_star+"/"+Run+"_round3_Unmapped.out.mate2"
    cmd+=" > "+options.output_star+"/"+Run+"_round4.output"
    cmd+=" 2> "+options.output_star+"/"+Run+"_round4.error"
    os.system(cmd)
    
    if options.star_shared_mem == True:
        # Remove a pre-loaded genome
        cmd  = "STAR "
        cmd += f" --runThreadN {options.cpu} "
        cmd += f" --genomeLoad Remove "
        cmd += f" --genomeDir {options.genome_dir_star} "
        os.system(cmd)
    
    cmd="rm "
    cmd+=options.output_star+"/"+Run+"_round4_Log.progress.out "+options.output_star+"/"+Run+"_round4_Log.out"
    os.system(cmd)  
    
    with logging_mutex:
        logger_proxy.info("STAR Round4 run for "+Run+" completed")     

def alignReadsWithOLegoRound5(options,Run,ended,condition,logger_proxy,logging_mutex,min_micro_exon): 
    if (os.path.exists(options.output_star+"/"+Run+"_olego_round5.sorted.bam")==True and samtoolsQuickCheck(options.output_star+"/"+Run+"_olego_round5.sorted.bam",options)==0) or os.path.exists(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam")==True:return 1
    cmd=options.softwares["olego"]
    cmd+=" -t "+str(options.cpu)
    cmd+=" -e "+str(min_micro_exon) # Minimum micro-exon size to be searched
    cmd+=" -M 0 " # No mismatches allowed
    cmd+=" -I 10000 " # Max intron size
    cmd+=" -i 20 " # Min intron size
    cmd+=" --max-multi 500 "
    cmd+=options.genome_dir_olego
    if ended=="SE":
        cmd+=" "+options.output_star+"/"+Run+"_round4_Unmapped.out.mate1"
        cmd+=" > "+options.output_star+"/"+Run+"_olego_round5.sam"
        cmd+=" 2> "+options.output_star+"/"+Run+"_round5.error"
        if os.path.exists(options.output_star+"/"+Run+"_olego_round5.sam")==False:
            os.system(cmd)
            
        cmd=options.softwares["xa2multi"]
        cmd+=" <("
        cmd+="cat "+options.output_star+"/"+Run+"_olego_round5.sam "
        cmd+="|awk '$2!=4') "
        cmd+=" > "+options.output_star+"/"+Run+"_olego_round5.modified.sam "
        #print(cmd)
        runCommand(["dummy",cmd])
        
        cmd="mv "
        cmd+=options.output_star+"/"+Run+"_olego_round5.modified.sam "
        cmd+=options.output_star+"/"+Run+"_olego_round5.sam "
        os.system(cmd)
        
        cmd="samtools "+" view "
        cmd+=" -@ "+str(options.cpu)
        cmd+=" -Sbh "
        cmd+=options.output_star+"/"+Run+"_olego_round5.sam "
        cmd+="|"
        cmd+="samtools "+" sort "
        cmd+=" -@ "+str(options.cpu)
        cmd+=" > "+options.output_star+"/"+Run+"_olego_round5.sorted.bam"
        os.system(cmd)
    else:
        cmd_f=cmd+" "+options.output_star+"/"+Run+"_round4_Unmapped.out.mate1"
        cmd_f+=" > "+options.output_star+"/"+Run+"_olego_round5_f.sam"
        cmd_f+=" 2> "+options.output_star+"/"+Run+"_round5_f.error"
        #print(cmd_f)
        if os.path.exists(options.output_star+"/"+Run+"_olego_round5_f.sam")==False:   
            os.system(cmd_f)
        
        cmd_r=cmd+" "+options.output_star+"/"+Run+"_round4_Unmapped.out.mate2"
        cmd_r+=" > "+options.output_star+"/"+Run+"_olego_round5_r.sam"
        cmd_r+=" 2> "+options.output_star+"/"+Run+"_round5_r.error"
        #print(cmd_r)
        if os.path.exists(options.output_star+"/"+Run+"_olego_round5_r.sam")==False:   
            os.system(cmd_r)
        
        cmd=options.softwares["mergePEsam.pl"]+" "
        cmd+=options.output_star+"/"+Run+"_olego_round5_f.sam "
        cmd+=options.output_star+"/"+Run+"_olego_round5_r.sam "
        cmd+=options.output_star+"/"+Run+"_olego_round5.sam"
        #print(cmd)
        if os.path.exists(options.output_star+"/"+Run+"_olego_round5.sam")==False:
            os.system(cmd)
        
        cmd=options.softwares["xa2multi"]
        cmd+=" <("
        cmd+="cat "+options.output_star+"/"+Run+"_olego_round5.sam "
        cmd+="|awk '$2!=4') "
        cmd+=" > "+options.output_star+"/"+Run+"_olego_round5.modified.sam "
        #print(cmd)
        runCommand(["dummy",cmd])
        
        cmd="mv "
        cmd+=options.output_star+"/"+Run+"_olego_round5.modified.sam "
        cmd+=options.output_star+"/"+Run+"_olego_round5.sam "
        os.system(cmd)
        
        cmd="samtools "+" view "
        cmd+=" -@ "+str(options.cpu)
        cmd+=" -Sbh "
        cmd+=options.output_star+"/"+Run+"_olego_round5.sam "
        cmd+="|"
        cmd+="samtools "+" sort "
        cmd+=" -@ "+str(options.cpu)
        cmd+=" > "+options.output_star+"/"+Run+"_olego_round5.sorted.bam"
        os.system(cmd)
    
    # Remove intermediate sam file
    os.system("rm -f "+options.output_star+"/"+Run+"_olego_round5.sam")
    os.system("rm -f "+options.output_star+"/"+Run+"_olego_round5_f.sam")
    os.system("rm -f "+options.output_star+"/"+Run+"_olego_round5_r.sam")
    with logging_mutex:
        logger_proxy.info("OLego run for "+Run+" completed")
    return 0
