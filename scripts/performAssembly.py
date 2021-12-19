##############################################################################################################################################################################
#
# List of functions for performing transcriptome assembly using PsiCLASS
#
##############################################################################################################################################################################

from scripts.alignReads import *
from scripts.fileReadWriteOperations import *
from scripts.runCommand import *
import multiprocessing

import pandas as pd


def reArrangeDataForAssembly( options, logger_proxy, logging_mutex ):
    all_runs = [os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/psiclass_output_sample_" + Run + ".gtf" ) for condition in options.mrna_md for Run in options.mrna_md[condition]]
    all_runs.append( os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined.gtf" ) )
    if False not in all_runs:return

    #########################################################################################################
    # Create list of introns
    #########################################################################################################
    pool = multiprocessing.Pool( processes = int( options.cpu ) )
    allinputs = []
    for condition_num, condition in enumerate( options.mrna_md ):
        for Run_num, Run in enumerate( options.mrna_md[condition] ):
            bam_filename_final = options.output_star + "/" + Run + "_final.sortedByCoord.out.bam"
            cmd = options.softwares["junc"]
            cmd += " " + bam_filename_final
            cmd += " -a "
            cmd += " > " + options.output_star + "/" + Run + "_introns"
            # os.system(cmd)
            if os.path.exists( options.output_star + "/" + Run + "_introns" ) == True:continue
            allinputs.append( [Run, cmd] )
    pool.map( runCommand, allinputs )

    for condition_num, condition in enumerate( options.mrna_md ):
        for Run_num, Run in enumerate( options.mrna_md[condition] ):
            cmd = "cat " + options.output_star + "/" + Run + "_introns|awk -vOFS=\"\\t\" '{print $1,$2+1,$3-1}' "
            cmd += " > " + options.output_star + "/" + Run + "_introns.bed"
            os.system( cmd )

    #########################################################################################################
    # Create list of exons
    #########################################################################################################
    pool = multiprocessing.Pool( processes = int( options.cpu ) )
    allinputs = []
    for condition_num, condition in enumerate( options.mrna_md ):
        for Run_num, Run in enumerate( options.mrna_md[condition] ):
            bam_filename_final = options.output_star + "/" + Run + "_final.sortedByCoord.out.bam"
            cmd = options.softwares["subexon-info"]
            cmd += " " + bam_filename_final
            cmd += " " + options.output_star + "/" + Run + "_introns"
            cmd += " --noStats "
            cmd += " > " + options.output_star + "/" + Run + "_exons"
            if os.path.exists( options.output_star + "/" + Run + "_exons" ) == True:continue
            allinputs.append( [Run, cmd] )
    pool.map( runCommand, allinputs )

    for condition_num, condition in enumerate( options.mrna_md ):
        for Run_num, Run in enumerate( options.mrna_md[condition] ):
            cmd = "cat " + options.output_star + "/" + Run + "_exons|grep -v ^#|awk -vOFS=\"\\t\" '{print $1,$2,$3}' "
            cmd += " > " + options.output_star + "/" + Run + "_exons.bed"
            if os.path.exists( options.output_star + "/" + Run + "_exons.bed" ) == True:continue
            os.system( cmd )

    pool = multiprocessing.Pool( processes = int( options.cpu ) )
    allinputs = []
    for condition_num, condition in enumerate( options.mrna_md ):
        for Run_num, Run in enumerate( options.mrna_md[condition] ):
            if os.path.exists( options.output_star + "/" + Run + "_num_exons_in_intron" ) == True:continue
            allinputs.append( [options.output_star + "/" + Run + "_exons.bed", options.output_star + "/" + Run + "_introns.bed", 0, options.output_star + "/" + Run + "_num_exons_in_intron"] )
    pool.map( findNumberOfExonsInEachIntron, allinputs )

    #########################################################################################################
    # Find the maximum number of exons housed within all introns less than 10K
    # Might be later substituted by a number
    #########################################################################################################
    threshold = 0
    for condition_num, condition in enumerate( options.mrna_md ):
        for Run_num, Run in enumerate( options.mrna_md[condition] ):
            filename = options.output_star + "/" + Run + "_num_exons_in_intron"
            fhr = open( filename, "r" )
            for line in fhr:
                chromosome, start, end, length, num_exons_in_introns = line.strip().split()
                length, num_exons_in_introns = int( length ), int( num_exons_in_introns )
                if threshold < num_exons_in_introns and length <= 10000:
                    threshold = num_exons_in_introns
            fhr.close()

    #########################################################################################################
    # Select the introns to be removed
    #########################################################################################################
    pool = multiprocessing.Pool( processes = int( options.cpu ) )
    allinputs = []
    for condition_num, condition in enumerate( options.mrna_md ):
        for Run_num, Run in enumerate( options.mrna_md[condition] ):
            filename = options.output_star + "/" + Run + "_num_exons_in_intron"
            prev_bamfilename = options.output_star + "/" + Run + "_final.sortedByCoord.out.bam"
            new_bamfilename = options.output_star + "/" + Run + "_for_psiclass.bam"
            if os.path.exists( new_bamfilename ) == True:continue
            allinputs.append( [filename,
                                       prev_bamfilename,
                                       new_bamfilename, threshold, options
                                       ] )
    pool.map( removeSpuriousLongIntronsInParallel, allinputs )


def runPsiCLASSMaxTerminalExonLength( options, logging_mutex, logger_proxy ):
    os.system( "mkdir -p " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined" )
    fhw = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/fofn", "w" )
    fhw_bg = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/bamgroup", "w" )
    run_to_counter = {}
    counter = 0
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            # fhw.write(options.output_star+"/"+Run+"_final.sortedByCoord.out.bam"+"\n")
            fhw.write( options.output_star + "/" + Run + "_for_psiclass.bam" + "\n" )
            fhw_bg.write( condition + "\n" )
            """if Run in run_to_counter:
                print("Repeated "+Run)"""
            run_to_counter[Run] = str( counter )
            # print(counter)
            counter += 1

    fhw.close()
    fhw_bg.close()
    os.system( "mkdir -p " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined" )
    # pprint.pprint(run_to_counter)
    # print(len(run_to_counter),counter)
    cmd = options.softwares["psiclass"] + " --lb "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/fofn "
    cmd += " -o " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/psiclass_output "
    cmd += " -p " + str( options.cpu )
    cmd += " --mateIdx 0 "
    cmd += " --bamGroup " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/bamgroup"
    cmd += " --primaryParalog "  # Outputs transcript assemblies with multi mapped reads
    cmd += " --tssTesQuantile 1"
    cmd += " > " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined.output "
    cmd += " 2> " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined.error "
    os.system( cmd )

    # Rename output files with Sample names
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            cmd = "mv " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/psiclass_output_sample_" + run_to_counter[Run] + ".gtf "
            cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/psiclass_output_sample_" + Run + ".gtf "
            os.system( cmd )
    cmd = "mv " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/psiclass_output_vote.gtf "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined.gtf "
    os.system( cmd )


def findNumberOfExonsInEachIntron( eachinput ):
    exon_bedfilename, intron_bedfilename, intron_min_length, outputfilename = eachinput
    exons = {}
    fhr = open( exon_bedfilename, "r" )
    for line in fhr:
        chromosome, start, end = line.strip().split()
        if chromosome not in exons:
            exons[chromosome] = []
        exons[chromosome].append( [int( start ), int( end )] )
    fhr.close()
    # pprint.pprint(exons)
    # sys.stdout.flush()
    exons_pd = {}
    for chromosome in exons:
        exons_pd[chromosome] = pd.DataFrame.from_dict( exons[chromosome] )
        exons_pd[chromosome].columns = ["start", "end"]
        # print(exons_pd[chromosome])
        # sys.stdout.flush()

    fhr = open( intron_bedfilename, "r" )
    fhw = open( outputfilename, "w" )
    for line in fhr:
        chromosome, start, end = line.strip().split()
        start, end = int( start ), int( end )
        if end - start + 1 < intron_min_length:continue
        if chromosome not in exons_pd:continue
        num_of_rows = exons_pd[chromosome][( exons_pd[chromosome].start >= start ) & ( exons_pd[chromosome].end <= end )].shape[0]
        # print(chromosome,start,end,end-start+1,num_of_rows)
        fhw.write( "\t".join( list( map( str, [chromosome, start, end, end - start + 1, num_of_rows] ) ) ) + "\n" )
        sys.stdout.flush()
    fhw.close()
    fhr.close()


def removeSpuriousLongIntronsInParallel( eachinput ):
    filename, prev_bamfilename, new_bamfilename, threshold, options = eachinput
    remove_these_introns = {}
    fhr = open( filename, "r" )
    for line in fhr:
        chromosome, start, end, length, num_exons_in_introns = line.strip().split()
        length, num_exons_in_introns = int( length ), int( num_exons_in_introns )
        if chromosome not in remove_these_introns:
            remove_these_introns[chromosome] = []
        if num_exons_in_introns > threshold:
            remove_these_introns[chromosome].append( start + "-" + end )
    fhr.close()

    """print(threshold)
    pprint.pprint(remove_these_introns)
    sys.stdout.flush()
    """
    cmd = "samtools view -h " + prev_bamfilename + " > " + prev_bamfilename[:-3] + "sam"
    if os.path.exists( prev_bamfilename[:-3] + "sam" ) == False:
        os.system( cmd )

    fhr = open( prev_bamfilename[:-3] + "sam", "r" )
    fhw = open( new_bamfilename[:-3] + "sam", "w" )
    for line in fhr:
        if line[0] == "@":
            fhw.write( line )
        else:
            chromosome = line.strip().split()[2]
            if chromosome not in remove_these_introns:
                fhw.write( line )
            elif "N" not in line.strip().split( "\t" )[5]:
                fhw.write( line )
            elif "RG:Z:4" not in line:
                fhw.write( line )
            else:
                coordinates = line.strip().split( "jI:B:i," )[-1].split()[0].split( "," )
                i = 0
                retain_this_alignment = 1
                while i < len( coordinates ):
                    junction_start, junction_end = coordinates[i], coordinates[i + 1]
                    if junction_start + "-" + junction_end in remove_these_introns[chromosome]:
                        retain_this_alignment = 0
                        break
                    i += 2
                if retain_this_alignment == 1:
                    fhw.write( line )
    fhw.close()
    fhr.close()
    cmd = "samtools view -Sbh " + new_bamfilename[:-3] + "sam > " + new_bamfilename
    os.system( cmd )

    """cmd="rm "+new_bamfilename[:-3]+"sam "+prev_bamfilename[:-3]+"sam"
    os.system(cmd)"""
