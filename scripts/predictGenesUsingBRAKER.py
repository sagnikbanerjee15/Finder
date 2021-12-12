
from pathlib import Path
from scripts.fileReadWriteOperations import *
from scripts.runCommand import *
import glob
import multiprocessing
import os
import time


def mapProteinsToGenomeUsingExonerate( protein_file, options, logger_proxy, logging_mutex ):
    """
    """
    all_proteins = readFastaFile( protein_file )
    # Split proteins into small files
    number_of_files = splitFasta( protein_file, protein_file[:-4], int( options.cpu ) - 1 )

    pool = multiprocessing.Pool( processes = int( options.cpu ) )
    allinputs = []
    for i in range( number_of_files ):
        cmd = "exonerate --model protein2genome --percent 90 --showcigar no -D 1000 "
        cmd += " --showquerygff  no --showtargetgff yes --showalignment no --showvulgar no --softmasktarget yes --softmaskquery yes "
        cmd += " --minintron 20 "
        # cmd+=" -c "+options.cpu+" "
        cmd += " -q " + f"{protein_file[:-4]}_{i}.fasta"
        cmd += " -t " + options.genome
        cmd += " > " + f"{options.output_assemblies_psiclass_terminal_exon_length_modified}/proteins_for_alignment_{i}.gff3 "
        cmd += " 2> " + f"{options.output_assemblies_psiclass_terminal_exon_length_modified}/proteins_for_alignment_{i}.error "
        with logging_mutex:
            logger_proxy.info( f"Running command - {cmd}" )
        allinputs.append( ["dummy", cmd] )
    pool.map( runCommand, allinputs )

    fhw = open( f"{options.output_assemblies_psiclass_terminal_exon_length_modified}/proteins_for_alignment.gff3", "w" )
    for i in range( number_of_files ):
        mapped_input_filename = f"{options.output_assemblies_psiclass_terminal_exon_length_modified}/proteins_for_alignment_{i}.gff3"
        fhr = open( mapped_input_filename, "r" )
        fhw.write( fhr.read() )
        fhr.close()
        os.system( f"rm {mapped_input_filename} {protein_file[:-4]}_{i}.fasta" )
    fhw.close()


def addBRAKERPredictions( options, logger_proxy, logging_mutex ):
    """
    """

    braker_gtf_file = options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.gtf"
    psiclass_gtf_file = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.gtf"
    if os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/FINDER_BRAKER_PROT.gtf" ) == True and os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_BRAKER_appended_high_conf.gtf" ) == True:return
    braker_transcripts = readAllTranscriptsFromGTFFileInParallel( [braker_gtf_file, "dummy", "dummy"] )[0]
    psiclass_transcripts = readAllTranscriptsFromGTFFileInParallel( [psiclass_gtf_file, "dummy", "dummy"] )[0]
    cmd = "perl -pe '/^>/ ? print \"\\n\" : chomp' "
    cmd += options.protein + "|tail -n +2 "
    cmd += "> " + options.protein + ".temp"
    os.system( cmd )

    cmd = "mv " + options.protein + ".temp " + options.protein
    os.system( cmd )
    all_proteins = readFastaFile( options.protein )

    cmd = "gffcompare "
    cmd += " -r " + braker_gtf_file
    cmd += " -o " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/braker_psiclass_gffcompare"
    cmd += " " + psiclass_gtf_file
    os.system( cmd )

    # remove_these=[]
    remove_these_braker_genes = []
    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/braker_psiclass_gffcompare.combined_with_CDS.gtf.tmap", "r" )
    for line in fhr:
        ref, transcript_id, ccode = line.strip().split()[:3]
        if transcript_id == "ref_id" or transcript_id == "-":continue
        # remove_these.append(transcript_id)
        remove_these_braker_genes.append( transcript_id.split( "." )[0] )
    fhr.close()
    # remove_these=list(set(remove_these))
    remove_these_braker_genes = list( set( remove_these_braker_genes ) )

    cmd = "makeblastdb "
    cmd += " -dbtype prot "
    cmd += " -in " + options.protein
    cmd += " -out " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/protein_db"
    os.system( cmd )

    cmd = "gffread "
    cmd += " " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.gtf "
    cmd += " -g " + options.genome
    cmd += " -w " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.fa "
    os.system( cmd )

    cmd = "perl -pe '/^>/ ? print \"\\n\" : chomp' "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.fa "
    cmd += "| tail -n +2 "
    cmd += " > " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.fa.temp "
    os.system( cmd )

    cmd = "mv " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.fa.temp "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.fa "
    os.system( cmd )

    cmd = "gffread "
    cmd += " " + psiclass_gtf_file
    cmd += " -g " + options.genome
    cmd += " -w " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.fasta"
    os.system( cmd )

    cmd = "perl -pe '/^>/ ? print \"\\n\" : chomp' "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.fasta"
    cmd += "| tail -n +2 "
    cmd += " > " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.fasta.temp"
    os.system( cmd )

    cmd = "mv " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.fasta.temp "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.fasta "
    os.system( cmd )

    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.fasta", "r" )
    fhw = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_prot.fasta", "w" )
    for line in fhr:
        if ">" in line:
            if "CDS" not in line:continue
            header, CDS = line.strip()[1:].split()
            CDS_start, CDS_end = list( map( int, CDS.split( "=" )[-1].split( "-" ) ) )
            seq = fhr.readline().strip()
            fhw.write( ">" + header + "\n" + translate( seq[CDS_start - 1:CDS_end] ) + "\n" )
    fhw.close()
    fhr.close()

    cmd = "blastp "
    cmd += " -db " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/protein_db"
    cmd += " -outfmt \"7 qseqid sseqid pident qcovs qcovhsp bitscore score\" "
    cmd += " -num_threads " + ( "32" if int( options.cpu ) >= 32 else options.cpu )
    cmd += " -out " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/finder_to_protein.out "
    cmd += " -query " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_prot.fasta "
    os.system( cmd )

    braker_fasta = readFastaFile( options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.fa" )
    all_braker_ids = list( braker_fasta.keys() )
    remove_these = []
    for each_braker_id in all_braker_ids:
        if each_braker_id.split( "." )[0] in remove_these_braker_genes:
            remove_these.append( each_braker_id )

    for id in remove_these:
        del braker_fasta[id]

    braker_proteins_fasta = {}
    for id in braker_fasta:
        braker_proteins_fasta[id] = translate( braker_fasta[id] )
    fhw = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker_prot.fa", "w" )
    for id in braker_proteins_fasta:
        fhw.write( ">" + id + "\n" + braker_proteins_fasta[id] + "\n" )

    fhw.close()

    cmd = "blastp "
    cmd += " -db " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/protein_db"
    cmd += " -outfmt \"7 qseqid sseqid pident qcovs qcovhsp bitscore score\" "
    cmd += " -num_threads " + ( "32" if int( options.cpu ) >= 32 else options.cpu )
    cmd += " -out " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker_to_protein.out "
    cmd += " -query " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker_prot.fa "
    os.system( cmd )

    cmd = "cat " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker_to_protein.out "
    cmd += "| grep -A 1 \"hits found\"|grep -v \"#\"|grep -v \"\\-\\-\"|awk '$3>=90 && $4>=90' "
    cmd += "|cut -f1 "
    cmd += "> "
    cmd += " " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker_to_protein_100_100.out "
    os.system( cmd )

    include_these_braker_models = []
    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker_to_protein_100_100.out", "r" )
    for line in fhr:
        include_these_braker_models.append( line.strip() )
    fhr.close()
    include_these_braker_models = set( include_these_braker_models )

    psiclass_and_braker_transcripts = readAllTranscriptsFromGTFFileInParallel( [options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_high_conf.gtf", "dummy", "dummy"] )[0]
    for chromosome in braker_transcripts:
        for transcript_id in braker_transcripts[chromosome]:
            if transcript_id in include_these_braker_models:
                if chromosome not in psiclass_and_braker_transcripts:
                    psiclass_and_braker_transcripts[chromosome] = {}
                psiclass_and_braker_transcripts[chromosome][transcript_id] = braker_transcripts[chromosome][transcript_id]

    # Find proteins that are not captured by either psiclass or braker
    protein_to_identity_coverage = {}
    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker_to_protein.out", "r" )
    for line in fhr:
        if line[0] == "#":continue
        gene, protein, identity, coverage, cov_per_hsp, bit_score, score = line.strip().split()
        identity = float( identity )
        coverage = float( coverage )
        if protein not in protein_to_identity_coverage:
            protein_to_identity_coverage[protein] = 0
        hm = 2 * identity * coverage / ( identity + coverage )
        if protein_to_identity_coverage[protein] <= hm:
            protein_to_identity_coverage[protein] = hm
    fhr.close()

    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/finder_to_protein.out", "r" )
    for line in fhr:
        if line[0] == "#":continue
        gene, protein, identity, coverage, cov_per_hsp, bit_score, score = line.strip().split()
        identity = float( identity )
        coverage = float( coverage )
        if protein not in protein_to_identity_coverage:
            protein_to_identity_coverage[protein] = 0
        hm = 2 * identity * coverage / ( identity + coverage )
        if protein_to_identity_coverage[protein] <= hm:
            protein_to_identity_coverage[protein] = hm
    fhr.close()
    proteins_for_alignment = {}
    for protein in protein_to_identity_coverage:
        # print(protein,protein_to_identity_coverage[protein])
        if protein_to_identity_coverage[protein] < 95:
            proteins_for_alignment[protein] = all_proteins[protein]

    fhw = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.fasta", "w" )
    for protein in proteins_for_alignment:
        fhw.write( ">" + protein + "\n" + proteins_for_alignment[protein] + "\n" )
    fhw.close()

    if options.exonerate_gff3 is None:
        mapProteinsToGenomeUsingExonerate( options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.fasta", options, logger_proxy, logging_mutex )
    else:
        cmd = f"cp {options.exonerate_gff3} {options.output_assemblies_psiclass_terminal_exon_length_modified}/proteins_for_alignment.gff3"
        os.system( cmd )

    cmd = options.softwares["convert_exonerate_gff_to_gtf"]
    cmd += " -i " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gff3 "
    cmd += " -o " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf "
    with logging_mutex:
        logger_proxy.info( f"Running command - {cmd}" )
    os.system( cmd )

    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf", "r" )
    fhw = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf.temp", "w" )
    for line in fhr:
        if line[0] == "#":
            fhw.write( line )
        else:
            chromosome = line.strip().split( "\t" )[0]
            fhw.write( line.replace( "gene_id \"", "gene_id \"" + chromosome + "." ).replace( "transcript_id \"", "transcript_id \"" + chromosome + "." ) )
    fhr.close()
    fhw.close()

    cmd = "mv " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf.temp "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf "
    os.system( cmd )

    cmd = "awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf "
    cmd += " > " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf.temp "
    os.system( cmd )

    cmd = "mv " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf.temp "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf "
    os.system( cmd )

    outputfilename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_BRAKER_appended_high_conf.gtf"
    writeTranscriptsToFile( [psiclass_and_braker_transcripts, outputfilename, 0] )

    cmd = "gffcompare "
    cmd += " -r " + outputfilename
    cmd += " -o " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/proteins_comparison_gffcompare "
    cmd += " " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf "
    os.system( cmd )

    proteins_to_be_discarded = []
    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_comparison_gffcompare.proteins_for_alignment.gtf.refmap", "r" )
    for line in fhr:
        ref_gene_id, ref_id, class_code, qry_id_list = line.strip().split()
        if "=" in class_code:
            proteins_to_be_discarded.append( qry_id_list.split( "|" )[0] )
    fhr.close()

    cmd = "cp " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_BRAKER_appended_high_conf.gtf "
    cmd += " " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/FINDER_BRAKER_PROT.gtf "
    os.system( cmd )

    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/proteins_for_alignment.gtf", "r" )
    fhw = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/FINDER_BRAKER_PROT.gtf", "a" )
    for line in fhr:
        if line.strip().split()[8].split( "gene_id" )[-1].split( ";" )[0].strip().strip( "\"" ) not in set( proteins_to_be_discarded ):
            fhw.write( line )
    fhw.close()
    fhr.close()

    writeTranscriptsToFile( [readAllTranscriptsFromGTFFileInParallel( [options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/FINDER_BRAKER_PROT.gtf", "dummy", "dummy"] )[0], options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/FINDER_BRAKER_PROT.gtf.temp", 0] )
    cmd = "cat " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_BRAKER_appended_high_conf.gtf "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_low_conf.gtf > "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_high_and_low_confidence_merged.gtf"
    os.system( cmd )


def configureAndRunBRAKER( options, logger_proxy, logging_mutex ):
    """
    """
    if os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.gtf" ) == True and Path( options.output_assemblies_psiclass_terminal_exon_length_modified + "/braker.gtf" ).stat().st_size != 0:return
    fhr = open( options.genome, "r" )
    fhw = open( options.temp_dir + "/genome_for_braker.fa", "w" )
    for line in fhr:
        if line[0] == ">":
            fhw.write( line.strip().split()[0] + "\n" )
        else:
            fhw.write( line )
    fhw.close()
    fhr.close()

    os.system( "rm -rf " + options.output_braker )
    os.system( "mkdir -p " + options.output_braker )
    """
    ################################################################################
    # Copy Augustus file to output directory
    ################################################################################
    cmd="cp -r "+options.softwares["augustus_main_dir"]+" "
    cmd+=options.output_braker+"/"
    os.system(cmd)
    with logging_mutex:
        logger_proxy.info(f"Running command - {cmd}")
    """

    with logging_mutex:
        logger_proxy.info( "BRAKER run started" )
    ################################################################################
    # Command to run BRAKER
    ################################################################################
    cores_for_braker = '40' if int( options.cpu ) > 40 else str( options.cpu )  # BRAKER cannot run with more than 40 CPUs
    cmd = options.softwares["braker"]
    cmd += " --GENEMARK_PATH=" + options.softwares["GENEMARK_PATH"]
    cmd += " --GUSHR_PATH=" + options.softwares["GUSHR_PATH"]
    cmd += " --softmasking "
    # cmd+=" --hints="+options.output_assemblies_psiclass_terminal_exon_length_modified+"/hints.gff "
    cmd += " --bam="
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            cmd += options.output_star + "/" + Run + "_for_psiclass.bam,"
    cmd += " --genome=" + options.temp_dir + "/genome_for_braker.fa "
    cmd += " --cores=" + cores_for_braker
    cmd += " --workingdir=" + options.output_braker
    cmd += " --overwrite "
    cmd += " --gff3 "
    if options.addUTR == True:
        cmd += " --addUTR=on "
    cmd += " > " + options.output_braker + ".output"
    cmd += " 2> " + options.output_braker + ".error"
    with logging_mutex:
        logger_proxy.info( f"Running BRAKER2 - {cmd}" )
    os.system( cmd )
    time.sleep( 5 )

    ################################################################################
    # Convert to gtf
    ################################################################################
    cmd = "gffread "
    cmd += " -E " + options.output_braker + "/braker.gff3"
    cmd += " -T -o " + options.output_braker + "/braker.gtf"
    os.system( cmd )
    with logging_mutex:
        logger_proxy.info( f"Running command - {cmd}" )

    cmd = "cp " + options.output_braker + "/braker.gtf "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/"
    os.system( cmd )

    if options.no_cleanup == False:
        # options.space_saved+=sum(f.stat().st_size for f in Path(options.output_braker).glob('**/*') if f.is_file() )
        # os.system("rm -rf "+options.output_braker)
        pass

