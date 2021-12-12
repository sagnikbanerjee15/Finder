
from scripts.fileReadWriteOperations import *
from scripts.runCommand import *
import multiprocessing
import os


def generateGenomicAndTranscriptomicCounts( options, logger_proxy, logging_mutex ):
    if os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_split_transcripts_with_bad_SJ_redundancy_removed.gtf" ) == True:return
    gtffilename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_transcripts_connecting_two_transcripts.gtf"
    writeTranscriptsToFile( [readAllTranscriptsFromGTFFileInParallel( [gtffilename, "dummy", "dummy"] )[0], gtffilename[:-4] + "_appended_info.gtf", 1] )
    exons_overlapping_with_introns_bedfilename = gtffilename[:-4] + "_exons_overlapping_with_introns.bed"
    cmd = "cat " + gtffilename[:-4] + "_appended_info.gtf |grep exon > " + gtffilename[:-4] + "_only_exons.gtf"
    # print(cmd)
    os.system( cmd )
    cmd = "cat " + gtffilename[:-4] + "_appended_info.gtf |grep intron > " + gtffilename[:-4] + "_only_introns.gtf"
    # print(cmd)
    os.system( cmd )
    cmd = "gtf2bed < " + gtffilename[:-4] + "_only_introns.gtf > " + gtffilename[:-4] + "_only_introns.bed"
    os.system( cmd )
    cmd = "gtf2bed < " + gtffilename[:-4] + "_only_exons.gtf > " + gtffilename[:-4] + "_only_exons.bed"
    os.system( cmd )
    cmd = "bedops --intersect " + gtffilename[:-4] + "_only_exons.bed " + gtffilename[:-4] + "_only_introns.bed > " + gtffilename[:-4] + "_exon_portions_overlapping_with_introns.bed"
    os.system( cmd )
    cmd = "bedops --element-of 1 " + gtffilename[:-4] + "_only_exons.bed " + gtffilename[:-4] + "_exon_portions_overlapping_with_introns.bed > " + exons_overlapping_with_introns_bedfilename
    os.system( cmd )
    pool = multiprocessing.Pool( processes = int( options.cpu ) )
    allinputs = []
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            if os.path.exists( options.output_star + "/" + Run + "_counts_all_info.pkl" ) == True:continue
            cmd = options.softwares["transferGenomicNucleotideCountsToTranscriptome"]
            cmd += " -b " + options.output_star + "/" + Run + "_for_psiclass.bam "
            # cmd+=" -g "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/"+condition+"/"+condition+"_appended.gtf"
            cmd += " -g " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_transcripts_connecting_two_transcripts.gtf "
            cmd += " -p " + options.output_star + "/" + Run + "_counts"
            cmd += " -f "
            cmd += " > " + options.output_star + "/" + Run + "_counts.output "
            cmd += " 2> " + options.output_star + "/" + Run + "_counts.error "

            # print(cmd)
            # os.system(cmd)
            # sys.exit()
            with logging_mutex:
                logger_proxy.info( cmd )
            allinputs.append( [Run, cmd] )
    pool.map( runCommand, allinputs )
