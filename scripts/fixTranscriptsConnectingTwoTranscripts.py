
from scripts.fileReadWriteOperations import *
import os


def reverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join([complement[base] for base in seq])[::-1]    

def fixTranscriptsConnectingTwoTranscripts(options,logging_mutex,logger_proxy):
    
    #########################################################################################################
    # Read all the transcript info
    #########################################################################################################
    
    
    complete_transcript_info={}
    exons_to_transcripts={}
    input_gtf_filename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_redundant_transcripts_removed.gtf"
    output_gtf_filename = options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_transcripts_connecting_two_transcripts.gtf"
    if os.path.exists(output_gtf_filename)==True:return
    dir_for_denovo_assemblies=options.output_assemblies_psiclass_terminal_exon_length_modified+"/denovo_assemblies"
    
    cmd="cp "
    cmd+=input_gtf_filename+" "
    cmd+=output_gtf_filename
    os.system(cmd)
    return []
    
    fhr=open(input_gtf_filename,"r")
    for line in fhr:
        chromosome,psiclass,structure,start,end,useless1,direction,useless2,desc=line.strip().split("\t")
        exon=start+"-"+end
        for ele in desc.split(";"):
            if "gene_id" in ele:
                gene_id=ele.split()[-1].strip("\"")
            if "transcript_id" in ele:
                transcript_id=ele.split()[-1].strip("\"")
        if structure=="transcript":
            continue
        if chromosome not in complete_transcript_info:
            complete_transcript_info[chromosome]={}
        if chromosome not in exons_to_transcripts:
            exons_to_transcripts[chromosome]={}
            
        if transcript_id not in complete_transcript_info[chromosome]:
            complete_transcript_info[chromosome][transcript_id]={"exons":[],
                                                                 "introns":[],
                                                                 "direction":direction,
                                                                 }
        if exon not in exons_to_transcripts[chromosome]:
            exons_to_transcripts[chromosome][exon]=[]
        if transcript_id not in exons_to_transcripts[chromosome][exon]:
            exons_to_transcripts[chromosome][exon].append(transcript_id)
        
        complete_transcript_info[chromosome][transcript_id]["exons"].append([int(start),int(end)])    
        
        # Construct the introns for each transcript
        if len(complete_transcript_info[chromosome][transcript_id]["exons"])>1:
            exon1=int(complete_transcript_info[chromosome][transcript_id]["exons"][-2][1])
            exon2=int(complete_transcript_info[chromosome][transcript_id]["exons"][-1][0])
            complete_transcript_info[chromosome][transcript_id]["introns"].append([exon1+1,exon2-1])
    fhr.close()
    
    
    #############################################################################################################
    # Find transcripts with two exons where both exons overlap with transcripts running in a different direction
    #############################################################################################################
    select_these_new_assemblies=[]
    all_regions_of_interest=[]
    all_transcripts_of_importance=[]
    regions_to_prev_transcripts_and_new_transcripts={}
    for chromosome in complete_transcript_info:
        for connecting_transcript_id in complete_transcript_info[chromosome]:
            if len(complete_transcript_info[chromosome][connecting_transcript_id]["exons"])==2:
                direction_of_connecting_transcript=complete_transcript_info[chromosome][connecting_transcript_id]["direction"]
                correct_overlap_of_exons=[0,0]
                transcripts_of_importance=[connecting_transcript_id]
                for exon_num,exon in enumerate(complete_transcript_info[chromosome][connecting_transcript_id]["exons"]):
                    for transcript_id in exons_to_transcripts[chromosome][str(exon[0])+"-"+str(exon[1])]:
                        if connecting_transcript_id!=transcript_id and direction_of_connecting_transcript!=complete_transcript_info[chromosome][transcript_id]["direction"] and direction_of_connecting_transcript!=complete_transcript_info[chromosome][transcript_id]["direction"]!="." and direction_of_connecting_transcript!=".":
                            correct_overlap_of_exons[exon_num]=1
                            transcripts_of_importance.append(transcript_id)
                if len(list(set(correct_overlap_of_exons)))==1 and list(set(correct_overlap_of_exons))[0]==1:
                    #print(condition,connecting_transcript_id)
                    transcripts_of_importance=list(set(transcripts_of_importance))
                    all_transcripts_of_importance.extend(transcripts_of_importance)
                    region_of_interest_start,region_of_interest_end=complete_transcript_info[chromosome][transcripts_of_importance[0]]["exons"][0][0],complete_transcript_info[chromosome][transcripts_of_importance[0]]["exons"][-1][1]
                    for transcript_id in transcripts_of_importance[1:]:
                        if region_of_interest_start>complete_transcript_info[chromosome][transcript_id]["exons"][0][0]:
                            region_of_interest_start=complete_transcript_info[chromosome][transcript_id]["exons"][0][0]
                        
                        if region_of_interest_end<complete_transcript_info[chromosome][transcript_id]["exons"][-1][1]:
                            region_of_interest_end=complete_transcript_info[chromosome][transcript_id]["exons"][-1][1]
                        
                    #print(condition,chromosome+":"+str(region_of_interest_start)+"-"+str(region_of_interest_end))
                    all_regions_of_interest.append(chromosome+":"+str(region_of_interest_start)+"-"+str(region_of_interest_end))
                    regions_to_prev_transcripts_and_new_transcripts[chromosome+":"+str(region_of_interest_start)+"-"+str(region_of_interest_end)]={"prev_transcripts":list(set(transcripts_of_importance)),
                                                                                                                                                   "new_transcripts":[]}
    #print("\n".join(list(set(all_transcripts_of_importance))))
    
    os.system("rm -rf "+dir_for_denovo_assemblies)
    os.system("mkdir -p "+dir_for_denovo_assemblies)
    #return
    rnaspades_transcript_number=1
    gene_number_SDN=0
    
    """
    # Generate index files
    pool = multiprocessing.Pool(processes=int(options.cpu))
    allinputs=[]
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            cmd="samtools index "+options.output_star+"/"+Run+"_for_psiclass.bam"
            #os.system(cmd)
            allinputs.append([Run,cmd])
            cmd="samtools index -c "+options.output_star+"/"+Run+"_for_psiclass.bam"
            #os.system(cmd)
            allinputs.append([Run,cmd])
    pool.map(runCommand,allinputs)
    """       
     
    for region_number,region in enumerate(all_regions_of_interest):
        with logging_mutex:
            logger_proxy.info("Starting de-novo assembly of region "+region)
        region_chromosome=region.split(":")[0]
        region_of_interest_start,region_of_interest_end=region.split(":")[-1].split("-")
        for condition in options.mrna_md:
            for Run in options.mrna_md[condition]:
                cmd="samtools "+" view "
                cmd+=" "+options.output_star+"/"+Run+"_for_psiclass.bam "
                cmd+=region
                if options.mrna_md[condition][Run]['Ended']=='SE':
                    cmd+="| grep -v ^@ | awk '{print \">\"$1\"\\n\"$10}' "
                    cmd+=" >> "+dir_for_denovo_assemblies+"/for_denovo_assembly_"+str(region_number)+".fasta"
                    os.system(cmd)
                else:
                    #cmd+="|grep -v ^@ | awk 'NR%2==1 {print \">\"$1\"_1\\n\"$10\"\\n}' "
                    cmd+="|grep -v ^@ | awk 'NR%2==1 {print \">\"$1\"_1\"\"\\n\"$10}'"
                    cmd+=" >> "+dir_for_denovo_assemblies+"/for_denovo_assembly_"+str(region_number)+".fasta"
                    os.system(cmd)
                    
                    cmd="samtools "+" view "
                    cmd+=" "+options.output_star+"/"+Run+"_for_psiclass.bam "
                    cmd+=region
                    #cmd+="|grep -v ^@ | awk 'NR%2==0 {print \">\"$1\"_2\\n\"$10\"\\n}' "
                    cmd+="|grep -v ^@ | awk 'NR%2==0 {print \">\"$1\"_1\"\"\\n\"$10}'"
                    cmd+=" >> "+dir_for_denovo_assemblies+"/for_denovo_assembly_"+str(region_number)+".fasta"
                    os.system(cmd)
            
        ###################################################################################
        # Generate de-novo assemblies
        ###################################################################################
        """cpu_for_spades=options.cpu if int(options.cpu)<30 else 30
        cmd="spades.py "
        cmd+=" --rna "
        cmd+=" --only-assembler "
        cmd+=" -t "+str(cpu_for_spades)
        cmd+=" -k 75 "
        cmd+=" --fast "
        cmd+=" -s "+dir_for_denovo_assemblies+"/for_denovo_assembly_"+str(region_number)+".fasta"
        cmd+=" -o "+dir_for_denovo_assemblies+"/rnaspades_"+str(region_number)
        cmd+=" > "+dir_for_denovo_assemblies+"/rnaspades_"+str(region_number)+".output "
        cmd+=" 2> "+dir_for_denovo_assemblies+"/rnaspades_"+str(region_number)+".error "
        os.system(cmd)
        with logging_mutex:
            logger_proxy.info("Command "+cmd+" completed")
        """
        transcript_lengths=[]
        for transcript_id in regions_to_prev_transcripts_and_new_transcripts[region]["prev_transcripts"]:
            transcript_lengths.append(sum([exon[1]-exon[0]+1 for exon in complete_transcript_info[region.split(":")[0]][transcript_id]["exons"]]))
        #print(transcript_lengths,sorted(transcript_lengths))
        #print(list(sorted(transcript_lengths))[1])
        desired_length=list(sorted(transcript_lengths))[1]+list(sorted(transcript_lengths))[2]
                
        # Creating configuration file for Soap Denovo Trans run
        fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/denovo_assemblies/config.txt","w")
        fhw.write("""#maximal read length
max_rd_len=100
[LIB]
#maximal read length in this lib
rd_len_cutof=45
#average insert size
avg_ins=200
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location
map_len=32
#fastq file for read 1
f="""+dir_for_denovo_assemblies+"/for_denovo_assembly_"+str(region_number)+".fasta")
        fhw.close()
        
        cmd="SOAPdenovo-Trans-127mer all "
        cmd+=" -s "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/denovo_assemblies/config.txt"
        cmd+=" -K 75 "
        cmd+=" -p "+(str(options.cpu) if int(options.cpu)<30 else str(30))
        cmd+=" -L "+str(int(0.95*desired_length))
        cmd+=" -o "+dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)
        cmd+=" > "+dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)+".output "
        cmd+=" 2> "+dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)+".error "
        os.system(cmd)
        with logging_mutex:
            logger_proxy.info("Command "+cmd+" completed")
        
        cmd="mv "
        cmd+=dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)+".scafSeq "
        cmd+=dir_for_denovo_assemblies+"/keepthis.scafSeq "
        os.system(cmd)
        
        cmd="rm "
        cmd+=dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)+".* "
        os.system(cmd)
        
        cmd="mv "
        cmd+=dir_for_denovo_assemblies+"/keepthis.scafSeq "
        cmd+=dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)+".scafSeq "
        os.system(cmd)
    
        ###################################################################################
        # Align de novo assemblies using GMAP
        ###################################################################################
        cmd="gmap "
        cmd+=" -t "+options.cpu
        cmd+=" -D "+options.genome_dir_gmap
        cmd+=" -d gmap_index "
        cmd+=" -f samse "
        cmd+=" --read-group-id=5 "
        cmd+=" --max-intronlength-middle=10000 --max-intronlength-ends=10000 "
        #cmd+=" "+dir_for_denovo_assemblies+"/"+str(region_number)+"/transcripts.fasta "    
        cmd+=" "+dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)+".scafSeq "
        cmd+=" > "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.sam"
        cmd+=" 2> "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly1.error"
        os.system(cmd)

        ###################################################################################
        # Remove uni-exon transcripts
        # Remove reads with multi-mappings 
        # Reverse complement reads based on orientation of alignment
        ###################################################################################
        
        fhr=open(dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.sam","r")
        fhw=open(dir_for_denovo_assemblies+"/"+str(region_number)+"_transcripts_altered.fasta","w")
        for line in fhr:
            if line[0]=="@":continue
            read_id,orientation,chromosome,pos,mapq,cigar,useless1,useless2,useless3,read_seq,read_qual=line.strip().split("\t")[:11]
            if "N" not in cigar:continue
            orientation=int(orientation)
            if orientation!=0 and orientation!=16:continue
            if (orientation==0 and "XS:A:-" in line) or (orientation==16 and "XS:A:+" in line):
                read_seq=reverseComplement(read_seq)
            fhw.write(">"+read_id+"\n"+read_seq+"\n")
        fhw.close()
        fhr.close()
        
        cmd="gmap "
        cmd+=" -t "+options.cpu
        cmd+=" -D "+options.genome_dir_gmap
        cmd+=" -d gmap_index "
        cmd+=" -f gff3_gene "
        cmd+=" --read-group-id=6 "
        cmd+=" --max-intronlength-middle=10000 --max-intronlength-ends=10000 "
        cmd+=" "+dir_for_denovo_assemblies+"/"+str(region_number)+"_transcripts_altered.fasta "
        cmd+=" > "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gff3"
        cmd+=" 2> "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly2.error"
        os.system(cmd)
        
        cmd="gffread "
        cmd+=" "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gff3"
        cmd+=" -T -o "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf"
        os.system(cmd)
        
        cmd="awk '$3!=\"CDS\"' "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf"
        cmd+=" > "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf.temp"
        os.system(cmd)
        
        cmd="mv "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf.temp "
        cmd+=dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf"
        os.system(cmd)
            
        ###################################################################################
        # Rename transcripts
        ###################################################################################
        fhr=open(dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf","r")
        soapdenovo_assembly=readAllTranscriptsFromGTFFileInParallel([dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf","dummy","dummy"])[0]
        genes_to_transcripts={}
        remove_these_chromosomes=[]
        for chromosome in soapdenovo_assembly:
            if chromosome!=region.split(":")[0]:
                remove_these_chromosomes.append(chromosome)
                continue
            for transcript_id in soapdenovo_assembly[chromosome]:
                gene_id=soapdenovo_assembly[chromosome][transcript_id]["gene_id"]
                if gene_id not in genes_to_transcripts:
                    genes_to_transcripts[gene_id]=[]
                genes_to_transcripts[gene_id].append(transcript_id)
        for chromosome in remove_these_chromosomes:
            del soapdenovo_assembly[chromosome]
            
        new_transcripts=[]
        #fhw=open(dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)+"_assembly.gtf.temp","w")
        new_name_gene_to_transcripts={}
        
        chromosome=region.split(":")[0]
        for gene_id in genes_to_transcripts:
            new_gene_name=chromosome+".SDN_"+str(gene_number_SDN)
            gene_number_SDN+=1
            new_name_gene_to_transcripts[new_gene_name]=[]
            for transcript_num,transcript_id in enumerate(genes_to_transcripts[gene_id]):
                new_transcript_id=chromosome+".SDN_"+str(gene_number_SDN)+"."+str(transcript_num)
                start,end=region.split(":")[-1].split("-")
                start=int(start)
                end=int(end)
                if soapdenovo_assembly[chromosome][transcript_id]["exons"][-1][1]<start or soapdenovo_assembly[chromosome][transcript_id]["exons"][0][0]>end:
                    del soapdenovo_assembly[chromosome][transcript_id]
                    continue
                new_transcripts.append(new_transcript_id)
                soapdenovo_assembly[chromosome][new_transcript_id]=soapdenovo_assembly[chromosome][transcript_id]
                soapdenovo_assembly[chromosome][new_transcript_id]["gene_id"]=new_gene_name
                del soapdenovo_assembly[chromosome][transcript_id]
        writeTranscriptsToFile([soapdenovo_assembly,dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf.temp"])
        fhr.close()
        
        cmd="cp "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf.temp "
        cmd+=dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gtf"
        os.system(cmd)    
        os.system("rm "+dir_for_denovo_assemblies+"/"+str(region_number)+"_assembly.gff3")
        
        regions_to_prev_transcripts_and_new_transcripts[region_chromosome+":"+str(region_of_interest_start)+"-"+str(region_of_interest_end)]["new_transcripts"]=list(set(new_transcripts))
        #print(region_chromosome+":"+str(region_of_interest_start)+"-"+str(region_of_interest_end))
        #return
    exclude_these=[]
    #print(len(regions_to_prev_transcripts_and_new_transcripts))
    #print(regions_to_prev_transcripts_and_new_transcripts)
    if len(regions_to_prev_transcripts_and_new_transcripts)>1:
        cmd="cat "+input_gtf_filename+" "
        for region_number,region in enumerate(all_regions_of_interest):
            #print(region_number,len(regions_to_prev_transcripts_and_new_transcripts[region]["new_transcripts"]))
            #print("Length of new region",len(regions_to_prev_transcripts_and_new_transcripts[region]["new_transcripts"]))
            with logging_mutex:
                logger_proxy.info(region+" "+str(len(regions_to_prev_transcripts_and_new_transcripts[region]["new_transcripts"]))+" "+";".join(regions_to_prev_transcripts_and_new_transcripts[region]["new_transcripts"]))
            if len(regions_to_prev_transcripts_and_new_transcripts[region]["new_transcripts"])==1:
                cmd+=dir_for_denovo_assemblies+"/soapdenovo_"+str(region_number)+"_assembly.gtf "
                exclude_these.extend(regions_to_prev_transcripts_and_new_transcripts[region]["prev_transcripts"])
        """
        if len(exclude_these)>0:
            cmd+="|grep -v \""+"\|".join([ele.replace(".","\.") for ele in exclude_these])+"\""
        """
        cmd+=" > "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_transcripts_connecting_two_transcripts.gtf "
        #print(cmd)
        os.system(cmd)
    else:
        cmd="cp "
        cmd+=input_gtf_filename+" "
        cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_transcripts_connecting_two_transcripts.gtf "
        os.system(cmd)
    
    #pprint.pprint(regions_to_prev_transcripts_and_new_transcripts)
    return regions_to_prev_transcripts_and_new_transcripts

