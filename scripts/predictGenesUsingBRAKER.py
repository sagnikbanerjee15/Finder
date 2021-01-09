
from pathlib import Path
from scripts.fileReadWriteOperations import *
import os
import time


def addBRAKERPredictions(options,logger_proxy,logging_mutex):
    """
    """
    
    braker_gtf_file=options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.gtf"
    psiclass_gtf_file=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS.gtf"
    
    
    #########################################################################################################################################
    # Split transcripts if BRAKER predicted a CDS in the UTRs of FINDER predictions
    #########################################################################################################################################
    '''
    # Fix direction of uniexon transcripts
    fhr=open(psiclass_gtf_file,"r")
    gene_to_direction={}
    for line in fhr:
        line=line.strip().split("\t")
        if line[2]!="transcript":continue
        gene_id=line[8].split("gene_id")[-1].split(";")[0].strip().strip("\"")
        #print(gene_id)
        direction=line[6]
        if gene_id not in gene_to_direction and direction!=".":
            gene_to_direction[gene_id]=direction
    #pprint.pprint(gene_to_direction)
    fhr.close()
    
    fhr=open(psiclass_gtf_file,"r")
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_direction_fixed.gtf","w")
    for line in fhr:
        line=line.strip().split("\t")
        gene_id=line[8].split("gene_id")[-1].split(";")[0].strip().strip("\"")
        if gene_id in gene_to_direction:
            line[6]=gene_to_direction[gene_id]
        fhw.write("\t".join(line)+"\n")
    
    fhw.close()
    fhr.close()
    
    cmd="gt gtf_to_gff3 "
    cmd+=" <(cat "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_direction_fixed.gtf"
    cmd+="|awk -vOFS=\"\t\" -F\"\t\" '$8=\".\"' "
    cmd+=")"
    cmd+="|sed 's/FPKM/fpkm/g'|sed 's/TPM/tpm/g'|gt gff3 -sort yes -retainids yes -tidy yes|"
    cmd+=options.softwares["canon-gff3"]
    cmd+=" > "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_direction_fixed.gff3 2> /dev/null"
    subprocess.check_call(['bash', '-c', cmd])
    
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_direction_fixed.gff3","r")
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_direction_fixed_only_UTR.gtf","w")
    for line in fhr:
        if line[0]=="#":
            #fhw.write(line)
            continue
        if line.strip().split("\t")[2]=="gene" or line.strip().split("\t")[2]=="mRNA":
            #fhw.write(line)
            if line.strip().split("\t")[2]=="mRNA":
                mRNA_def=";".join(line.strip().split("\t")[-1].split(";")[2:])
                for ele in mRNA_def.split(";"):
                    if "gene_id" in ele:
                        gene_id=ele.split("gene_id=")[-1]
                    if "transcript_id" in ele:
                        transcript_id_num=ele.split("transcript_id=")[-1].split(".")[-1]
        if line.strip().split("\t")[2]=="five_prime_UTR":
            #line=line.strip()+";"+mRNA_def
            line=line.split("\t")
            line[2]="exon"
            line[-1]="transcript_id \""+gene_id+"_5_prime_UTR."+transcript_id_num+"\"; gene_id \""+gene_id+"_5_prime_UTR\""
            line="\t".join(line)
            
            fhw.write(line+"\n")
            
        if line.strip().split("\t")[2]=="three_prime_UTR":
            #line=line.strip()+";"+mRNA_def
            line=line.split("\t")
            line[2]="exon"
            line[-1]="transcript_id \""+gene_id+"_3_prime_UTR."+transcript_id_num+"\"; gene_id \""+gene_id+"_3_prime_UTR\""
            line="\t".join(line)
            fhw.write(line+"\n")
    fhw.close()
    fhr.close()
    
    cmd="mikado compare "
    cmd+=" -r "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_direction_fixed_only_UTR.gtf"
    cmd+=" -p "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.gtf"
    cmd+=" -o "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/braker_in_finder_UTRs_mikado"
    cmd+=" -erm "
    os.system(cmd)
    
    
    braker_transcripts=[]
    finder_transcripts=[]
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/braker_in_finder_UTRs_mikado.refmap","r")
    for line in fhr:
        finder_transcript,ccode,braker_transcript=line.strip().split("\t")[:3]
        if ccode=="=" or ccode=="_":
            braker_transcripts.append(braker_transcript)
            finder_transcripts.append(finder_transcript)
    fhr.close()
    
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.gtf","r")
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_modify_finder_whole_transcript.gtf","w")
    fhw_CDS=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_modify_finder_CDS.gtf","w")
    
    for line in fhr:
        if line.strip().split("\t")[2]=="transcript":continue
        if line.strip().split("transcript_id")[-1].split(";")[0].strip().strip("\"") in braker_transcripts:
            fhw.write(line)
            if line.strip().split("\t")[2]=="CDS":
                fhw_CDS.write(line)
    fhr.close()
    fhw.close()
    fhw_CDS.close()
    
    #print(finder_transcripts)
    finder_transcript_modified=[]
    for transcript in finder_transcripts:
        if "_5_prime_UTR" in transcript:
            finder_transcript_modified.append(transcript.replace("_5_prime_UTR",''))
            #print(transcript.replace("_5_prime_UTR",''))
        if "_3_prime_UTR" in transcript:
            finder_transcript_modified.append(transcript.replace("_3_prime_UTR",''))
            #print(transcript.replace("_3_prime_UTR",''))
    finder_transcripts=finder_transcript_modified
    #print(finder_transcripts)
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_direction_fixed.gtf","r")
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript.gtf","w")
    fhw_CDS=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined_with_CDS_direction_fixed_to_modify_finder_CDS.gtf","w")
    for line in fhr:
        if line.strip().split("\t")[2]=="transcript":continue
        if line.strip().split("transcript_id")[-1].split(";")[0].strip().strip("\"") in finder_transcripts:
            fhw.write(line)
            if line.strip().split("\t")[2]=="CDS":
                fhw_CDS.write(line)
        
    fhr.close()
    fhw.close()
    fhw_CDS.close()
    
    combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript=readAllTranscriptsFromGTFFileInParallel([options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript.gtf","dummy","dummy"])[0]
    temp={}
    for chromosome in combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript:
        for transcript_id in combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript[chromosome]:
            temp[transcript_id]=combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript[chromosome][transcript_id]
    combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript=temp
    
    fhw_first_and_last_exon=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/finder_first_and_last_exon_point.gtf","w")
    for transcript_id in combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript:
        """
        whole_annotations[chromosome][transcript_id]={"exons":[],
                                                        "introns":[],
                                                        "cds":[],
                                                        "cds_frame":[],
                                                        "direction":direction,
                                                        "TPM":tpm,
                                                        "cov":cov,
                                                        "gene_id":gene_id,
                                                        "FPKM":fpkm,
                                                        "transcript_start":0,
                                                        "transcript_end":0,
                                                        "chromosome":chromosome
                                                        }
        """
        
        chromosome=combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript[transcript_id]["chromosome"]
        starting_exon=combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript[transcript_id]["exons"][0]
        last_exon=combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript[transcript_id]["exons"][-1]
        gene_id=combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript[transcript_id]["gene_id"]
        direction=combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript[transcript_id]["direction"]
        
        line=[chromosome,"PsiCLASS","exon",str(starting_exon[0]),str(starting_exon[0]+1),"1000",direction,".","gene_id \""+gene_id+"\"; transcript_id \""+transcript_id+"\""]
        fhw_first_and_last_exon.write("\t".join(line)+"\n")
        
        line=[chromosome,"PsiCLASS","exon",str(last_exon[1]-1),str(last_exon[1]),"1000",direction,".","gene_id \""+gene_id+"\"; transcript_id \""+transcript_id+"\""]
        fhw_first_and_last_exon.write("\t".join(line)+"\n")
    fhw_first_and_last_exon.close()
    
    cmd="bedtools subtract -A "
    cmd+=" -a <(bedtools subtract -a "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined_with_CDS_direction_fixed_to_modify_finder_whole_transcript.gtf"+" -b <(cat "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_modify_finder_CDS.gtf "+ options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined_with_CDS_direction_fixed_to_modify_finder_CDS.gtf ) ) "
    cmd+=" -b "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/finder_first_and_last_exon_point.gtf "
    cmd+=" > "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/exons_to_be_targeted.gtf"
    subprocess.check_call(['bash', '-c', cmd])
    '''
    #########################################################################################################################################
    
    
    """            
    exons_to_be_targeted=readAllTranscriptsFromGTFFileInParallel([options.output_assemblies_psiclass_terminal_exon_length_modified+"/exons_to_be_targeted.gtf","dummy","dummy"])[0]
    FINDER_BRAKER_PROT=readAllTranscriptsFromGTFFileInParallel([psiclass_gtf_file,"dummy","dummy"])[0]
                                                                  
    #pprint.pprint(exons_to_be_targeted)
    sys.stdout.flush()
    add_these_transcripts={}
    remove_these_transcripts={}
    for chromosome in exons_to_be_targeted:
        for transcript_id in exons_to_be_targeted[chromosome]:
            if len(exons_to_be_targeted[chromosome][transcript_id]["exons"])==1:
                exon_start,exon_end=exons_to_be_targeted[chromosome][transcript_id]["exons"][0][0],exons_to_be_targeted[chromosome][transcript_id]["exons"][0][1]
                modified_exon_start,modified_exon_end=int((exon_start+exon_end)/2)-10,int((exon_start+exon_end)/2)+10
                print(transcript_id,chromosome+":"+str(int((exon_start+exon_end)/2-10))+"-"+str(int((exon_start+exon_end)/2+10)))
                for transcript_id_to_be_split in FINDER_BRAKER_PROT[chromosome]:
                    transcript_exon_start,transcript_exon_end=FINDER_BRAKER_PROT[chromosome][transcript_id_to_be_split]["transcript_start"],FINDER_BRAKER_PROT[chromosome][transcript_id_to_be_split]["transcript_end"]
                    if transcript_exon_start<=exon_start<=transcript_exon_end and transcript_exon_start<=exon_end<=transcript_exon_end:
                        
                        for exon in FINDER_BRAKER_PROT[chromosome][transcript_id_to_be_split]["exons"]:
                            if exon[0]<exon_start<exon[1] and exon[0]<exon_end<exon[1]:
                                #print(transcript_id_to_be_split,chromosome+":"+str(exon[0])+"-"+str(exon[1]))
                                FINDER_BRAKER_PROT[chromosome][transcript_id_to_be_split]
                                new_transcript_id_1=chromosome+"."+transcript_id_to_be_split.split(".")[1]+"_1b."+transcript_id_to_be_split.split(".")[-1]
                                new_transcript_id_2=chromosome+"."+transcript_id_to_be_split.split(".")[1]+"_2b."+transcript_id_to_be_split.split(".")[-1]
                                new_gene_id_1=".".join(new_transcript_id_1.split(".")[:2])
                                new_gene_id_2=".".join(new_transcript_id_2.split(".")[:2])
                                fill_transcript_id_1=1
                                direction=FINDER_BRAKER_PROT[chromosome][transcript_id_to_be_split]["direction"]
                                if chromosome not in add_these_transcripts:
                                    add_these_transcripts[chromosome]={}
                                add_these_transcripts[chromosome][new_transcript_id_1]={"exons":[],
                                                        "introns":[],
                                                        "cds":[],
                                                        "cds_frame":[],
                                                        "direction":direction,
                                                        "TPM":0,
                                                        "cov":0,
                                                        "gene_id":new_gene_id_1,
                                                        "FPKM":0,
                                                        "transcript_start":0,
                                                        "transcript_end":0,
                                                        "chromosome":chromosome}
                                for exon in FINDER_BRAKER_PROT[chromosome][transcript_id_to_be_split]["exons"]:
                                    if exon[0]<exon_start<exon[1] and exon[0]<exon_end<exon[1]:
                                        fill_transcript_id_1=0
                                        
                                        add_these_transcripts[chromosome][new_transcript_id_1]["exons"].append([exon[0],int((exon_start+exon_end)/2)-10])
                                        #add_these_transcripts[chromosome][new_transcript_id_1]["exons"][-1][1]=int((exon_start+exon_end)/2)-10
                                        add_these_transcripts[chromosome][new_transcript_id_2]={"exons":[],
                                                        "introns":[],
                                                        "cds":[],
                                                        "cds_frame":[],
                                                        "direction":direction,
                                                        "TPM":0,
                                                        "cov":0,
                                                        "gene_id":new_gene_id_2,
                                                        "FPKM":0,
                                                        "transcript_start":0,
                                                        "transcript_end":0,
                                                        "chromosome":chromosome}
                                        if transcript_id_to_be_split=="1.3193_opp.0":
                                            print("Adjusting",chromosome+":"+str(exon_start)+"-"+str(exon_end))
                                            print("New Exon",chromosome+":"+str(int((exon_start+exon_end)/2)+10)+"-"+str(exon[1]))
                                        
                                        if transcript_id_to_be_split=="1.3193_opp.0":
                                            print("Writing ",exon," to second split")
                                        add_these_transcripts[chromosome][new_transcript_id_2]["exons"].append([int((exon_start+exon_end)/2)+10,exon[1]])
                                        #add_these_transcripts[chromosome][new_transcript_id_2]["exons"][0][0]=int((exon_start+exon_end)/2)+10
                                        continue
                                    if fill_transcript_id_1==1:
                                        if transcript_id_to_be_split=="1.3193_opp.0":
                                            print("Writing ",exon," to first split")
                                        add_these_transcripts[chromosome][new_transcript_id_1]["exons"].append(exon)
                                    else:
                                        if transcript_id_to_be_split=="1.3193_opp.0":
                                            print("Writing ",exon," to second split")
                                        add_these_transcripts[chromosome][new_transcript_id_2]["exons"].append(exon)
                                if chromosome not in remove_these_transcripts:
                                    remove_these_transcripts[chromosome]=[]
                                remove_these_transcripts[chromosome].append(transcript_id_to_be_split)
                                break
                            
    for chromosome in remove_these_transcripts:
        for transcript in set(remove_these_transcripts[chromosome]):
            del FINDER_BRAKER_PROT[chromosome][transcript]
    
    for chromosome in add_these_transcripts:
        for transcript in add_these_transcripts[chromosome]:
            FINDER_BRAKER_PROT[chromosome][transcript]=add_these_transcripts[chromosome][transcript]
    writeTranscriptsToFile([FINDER_BRAKER_PROT,options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_models_fixed.gtf"])
    return
    """
    
    
    #return
    if os.path.exists(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/FINDER_BRAKER_PROT.gtf")==True and os.path.exists(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_BRAKER_appended_high_conf.gtf")==True:return
    braker_transcripts=readAllTranscriptsFromGTFFileInParallel([braker_gtf_file,"dummy","dummy"])[0]
    psiclass_transcripts=readAllTranscriptsFromGTFFileInParallel([psiclass_gtf_file,"dummy","dummy"])[0]
    cmd="perl -pe '/^>/ ? print \"\\n\" : chomp' "
    cmd+=options.protein+"|tail -n +2 "
    cmd+="> "+options.protein+".temp"
    os.system(cmd)
        
    cmd="mv "+options.protein+".temp "+options.protein
    os.system(cmd)
    all_proteins=readFastaFile(options.protein)
    
    cmd="gffcompare "
    cmd+=" -r "+braker_gtf_file
    cmd+=" -o "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/braker_psiclass_gffcompare"
    cmd+=" "+psiclass_gtf_file
    os.system(cmd)
    
    #remove_these=[]
    remove_these_braker_genes=[]
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/braker_psiclass_gffcompare.combined_with_CDS.gtf.tmap","r")
    for line in fhr:
        ref,transcript_id,ccode=line.strip().split()[:3]
        if transcript_id=="ref_id" or transcript_id=="-":continue
        #remove_these.append(transcript_id)
        remove_these_braker_genes.append(transcript_id.split(".")[0])
    fhr.close()
    #remove_these=list(set(remove_these))
    remove_these_braker_genes=list(set(remove_these_braker_genes))
    
    cmd="makeblastdb "
    cmd+=" -dbtype prot "
    cmd+=" -in "+options.protein
    cmd+=" -out "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/protein_db"
    os.system(cmd)
    
    cmd="gffread "
    cmd+=" "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.gtf "
    cmd+=" -g "+options.genome
    cmd+=" -w "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.fa "
    os.system(cmd)
    
    cmd="perl -pe '/^>/ ? print \"\\n\" : chomp' "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.fa "
    cmd+="| tail -n +2 "
    cmd+=" > "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.fa.temp "
    os.system(cmd)
    
    cmd="mv "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.fa.temp "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.fa "
    os.system(cmd)

    cmd="gffread "
    cmd+=" "+psiclass_gtf_file
    cmd+=" -g "+options.genome
    cmd+=" -w "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS.fasta"
    os.system(cmd)
    
    cmd="perl -pe '/^>/ ? print \"\\n\" : chomp' "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS.fasta"
    cmd+="| tail -n +2 "
    cmd+=" > "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS.fasta.temp"
    os.system(cmd)
    
    cmd="mv "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS.fasta.temp "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS.fasta "
    os.system(cmd)
    
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS.fasta","r")
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_prot.fasta","w")
    for line in fhr:
        if ">" in line:
            if "CDS" not in line:continue
            header,CDS=line.strip()[1:].split()
            CDS_start,CDS_end=list(map(int,CDS.split("=")[-1].split("-")))
            seq=fhr.readline().strip()
            fhw.write(">"+header+"\n"+translate(seq[CDS_start-1:CDS_end])+"\n")
    fhw.close()
    fhr.close()
    
    cmd="blastp "
    cmd+=" -db "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/protein_db"
    cmd+=" -outfmt \"7 qseqid sseqid pident qcovs qcovhsp bitscore score\" "
    cmd+=" -num_threads "+ ("32" if int(options.cpu)>=32 else options.cpu )
    cmd+=" -out "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/finder_to_protein.out "
    cmd+=" -query "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_prot.fasta "
    if os.path.exists(options.output_assemblies_psiclass_terminal_exon_length_modified+"/finder_to_protein.out")==False:
        os.system(cmd)
    
    braker_fasta=readFastaFile(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.fa")
    all_braker_ids=list(braker_fasta.keys())
    remove_these=[]
    for each_braker_id in all_braker_ids:
        if each_braker_id.split(".")[0] in remove_these_braker_genes:
            remove_these.append(each_braker_id)
    
    for id in remove_these:
        del braker_fasta[id]
    
    braker_proteins_fasta={}
    for id in braker_fasta:
        braker_proteins_fasta[id]=translate(braker_fasta[id])
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_prot.fa","w")
    for id in braker_proteins_fasta:
        fhw.write(">"+id+"\n"+braker_proteins_fasta[id]+"\n")
    
    fhw.close()
    
    cmd="blastp "
    cmd+=" -db "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/protein_db"
    cmd+=" -outfmt \"7 qseqid sseqid pident qcovs qcovhsp bitscore score\" "
    cmd+=" -num_threads "+ ("32" if int(options.cpu)>=32 else options.cpu )
    cmd+=" -out "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_protein.out "
    cmd+=" -query "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_prot.fa "
    if os.path.exists(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_protein.out")==False:
        os.system(cmd)

    cmd="cat "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_protein.out "
    cmd+="| grep -A 1 \"hits found\"|grep -v \"#\"|grep -v \"\\-\\-\"|awk '$3>=90 && $4>=90' "
    cmd+="|cut -f1 "
    cmd+="> "
    cmd+=" "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_protein_100_100.out "
    os.system(cmd)
    
    include_these_braker_models=[]
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_protein_100_100.out","r")
    for line in fhr:
        include_these_braker_models.append(line.strip())
    fhr.close()
    include_these_braker_models=set(include_these_braker_models)
    
    psiclass_and_braker_transcripts=readAllTranscriptsFromGTFFileInParallel([options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_high_conf.gtf","dummy","dummy"])[0]
    for chromosome in braker_transcripts:
        for transcript_id in braker_transcripts[chromosome]:
            if transcript_id in include_these_braker_models:
                if chromosome not in psiclass_and_braker_transcripts:
                    psiclass_and_braker_transcripts[chromosome]={}
                psiclass_and_braker_transcripts[chromosome][transcript_id]=braker_transcripts[chromosome][transcript_id]
    
    # Find proteins that are not captured by either psiclass or braker
    protein_to_identity_coverage={}
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker_to_protein.out","r")
    for line in fhr:
        if line[0]=="#":continue
        gene,protein,identity,coverage,cov_per_hsp,bit_score,score=line.strip().split()
        identity=float(identity)
        coverage=float(coverage)
        if protein not in protein_to_identity_coverage:
            protein_to_identity_coverage[protein]=0
        hm=2*identity*coverage/(identity+coverage)
        if protein_to_identity_coverage[protein]<=hm:
            protein_to_identity_coverage[protein]=hm
    fhr.close()
    
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/finder_to_protein.out","r")
    for line in fhr:
        if line[0]=="#":continue
        gene,protein,identity,coverage,cov_per_hsp,bit_score,score=line.strip().split()
        identity=float(identity)
        coverage=float(coverage)
        if protein not in protein_to_identity_coverage:
            protein_to_identity_coverage[protein]=0
        hm=2*identity*coverage/(identity+coverage)
        if protein_to_identity_coverage[protein]<=hm:
            protein_to_identity_coverage[protein]=hm
    fhr.close()
    proteins_for_alignment={}
    for protein in protein_to_identity_coverage:
        #print(protein,protein_to_identity_coverage[protein])
        if protein_to_identity_coverage[protein]<95:
            proteins_for_alignment[protein]=all_proteins[protein]
    
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.fasta","w")
    for protein in proteins_for_alignment:
        fhw.write(">"+protein+"\n"+proteins_for_alignment[protein]+"\n")
    fhw.close()
    
    cmd="exonerate --model protein2genome --percent 90 --showcigar no -D 1000 "
    #cmd+=" --exhaustive "
    cmd+=" --showquerygff  no --showtargetgff yes --showalignment no --showvulgar no --softmasktarget yes "
    cmd+=" --minintron 20 " 
    cmd+=" -c "+options.cpu+" "
    cmd+=" -q "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.fasta "
    cmd+=" -t "+options.genome
    cmd+=" > "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gff3 "
    os.system(cmd)
    
    cmd=options.softwares["convert_exonerate_gff_to_gtf"]
    cmd+=" -i "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gff3 "
    cmd+=" -o "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf "
    os.system(cmd)
    
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf","r")
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf.temp","w")
    for line in fhr:
        if line[0]=="#":
            fhw.write(line)
        else:
            chromosome=line.strip().split("\t")[0]
            fhw.write(line.replace("gene_id \"","gene_id \""+chromosome+".").replace("transcript_id \"","transcript_id \""+chromosome+"."))
    fhr.close()
    fhw.close()
    
    cmd="mv "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf.temp "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf "
    os.system(cmd)
     
    cmd="awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf "
    cmd+=" > "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf.temp "
    os.system(cmd)
    
    cmd="mv "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf.temp "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf "
    os.system(cmd)
    
    outputfilename=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_BRAKER_appended_high_conf.gtf"
    writeTranscriptsToFile([psiclass_and_braker_transcripts,outputfilename,0])

    cmd="gffcompare "
    cmd+=" -r "+outputfilename
    cmd+=" -o "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/proteins_comparison_gffcompare "
    cmd+=" "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf "
    os.system(cmd)
    
    proteins_to_be_discarded=[]
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_comparison_gffcompare.proteins_for_alignment.gtf.refmap","r")
    for line in fhr:
        ref_gene_id,ref_id,class_code,qry_id_list=line.strip().split()
        if "=" in class_code:
            proteins_to_be_discarded.append(qry_id_list.split("|")[0])
    fhr.close()
    
    cmd="cp "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_BRAKER_appended_high_conf.gtf "
    cmd+=" "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/FINDER_BRAKER_PROT.gtf "
    os.system(cmd)
    
    fhr=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/proteins_for_alignment.gtf","r")
    fhw=open(options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/FINDER_BRAKER_PROT.gtf","a")
    for line in fhr:
        if line.strip().split()[8].split("gene_id")[-1].split(";")[0].strip().strip("\"") not in set(proteins_to_be_discarded):
            fhw.write(line)
    fhw.close()
    fhr.close()
    
    writeTranscriptsToFile([readAllTranscriptsFromGTFFileInParallel([options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/FINDER_BRAKER_PROT.gtf","dummy","dummy"])[0],options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/FINDER_BRAKER_PROT.gtf.temp",0])
    cmd="cat "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_BRAKER_appended_high_conf.gtf "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_low_conf.gtf > "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_with_CDS_high_and_low_confidence_merged.gtf"
    os.system(cmd)
       
def configureAndRunBRAKER(options,logger_proxy,logging_mutex):
    """
    """
    if os.path.exists(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.gtf")==True and Path(options.output_assemblies_psiclass_terminal_exon_length_modified+"/braker.gtf").stat().st_size!=0:return
    fhr=open(options.genome,"r")
    fhw=open(options.temp_dir+"/genome_for_braker.fa","w")
    for line in fhr:
        if line[0]==">":
            fhw.write(line.strip().split()[0]+"\n")
        else:
            fhw.write(line)
    fhw.close()
    fhr.close()
    
    os.system("rm -rf "+options.output_braker)
    os.system("mkdir -p "+options.output_braker)
    
    ################################################################################
    # Copy Augustus file to output directory
    ################################################################################
    cmd="cp -r "+options.softwares["augustus_main_dir"]+" "
    cmd+=options.output_braker+"/"
    os.system(cmd)    
    with logging_mutex:
        logger_proxy.info(f"Running command - {cmd}")
    
    with logging_mutex:
        logger_proxy.info("BRAKER run started")
    ################################################################################
    # Command to run BRAKER
    ################################################################################
    cores_for_braker='40' if int(options.cpu)>40 else str(options.cpu) # BRAKER cannot run with more than 40 CPUs
    cmd=options.softwares["braker"]
    cmd+=" --GENEMARK_PATH="+options.softwares["GENEMARK_PATH"]
    cmd+=" --AUGUSTUS_CONFIG_PATH="+options.softwares["AUGUSTUS_CONFIG_PATH"]
    cmd+=" --AUGUSTUS_BIN_PATH="+options.softwares["AUGUSTUS_BIN_PATH"]
    cmd+=" --AUGUSTUS_SCRIPTS_PATH="+options.softwares["AUGUSTUS_SCRIPTS_PATH"]
    cmd+=" --GUSHR_PATH="+options.softwares["GUSHR_PATH"]
    cmd+=" --softmasking "
    #cmd+=" --hints="+options.output_assemblies_psiclass_terminal_exon_length_modified+"/hints.gff "
    cmd+=" --bam="
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            cmd+=options.output_star+"/"+Run+"_for_psiclass.bam,"
    cmd+=" --genome="+options.temp_dir+"/genome_for_braker.fa "
    cmd+=" --cores="+cores_for_braker
    cmd+=" --workingdir="+options.output_braker
    cmd+=" --overwrite "
    cmd+=" --gff3 "
    cmd+=" --addUTR=on "
    #cmd+=" --UTR=on "
    cmd+=" > "+options.output_braker+".output"
    cmd+=" 2> "+options.output_braker+".error"
    #print(cmd)
    with logging_mutex:
        logger_proxy.info(f"Running BRAKER2 - {cmd}")
    os.system(cmd)
    time.sleep(5)

    ################################################################################
    # Convert to gtf
    ################################################################################
    cmd="gffread "
    cmd+=" -E "+options.output_braker+"/braker.gff3"
    cmd+=" -T -o "+options.output_braker+"/braker.gtf"
    os.system(cmd)
    with logging_mutex:
        logger_proxy.info(f"Running command - {cmd}")
    
    cmd="cp "+options.output_braker+"/braker.gtf "
    cmd+=options.output_assemblies_psiclass_terminal_exon_length_modified+"/"
    os.system(cmd)
    
    if options.no_cleanup==False:
        #options.space_saved+=sum(f.stat().st_size for f in Path(options.output_braker).glob('**/*') if f.is_file() )
        #os.system("rm -rf "+options.output_braker)
        pass
 
