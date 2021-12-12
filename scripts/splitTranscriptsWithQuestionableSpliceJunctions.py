

from copy import deepcopy
from scripts.fileReadWriteOperations import *
from scripts.removeRedundantTranscripts import *
import multiprocessing

import numpy as np


def extractSJ( eachinput ):
    """
    """
    bam_filename, output_filename, options = eachinput
    cmd = "regtools "
    cmd += " junctions "
    cmd += " extract "
    cmd += " -a 4 "
    cmd += " -m 20 "
    cmd += " -M 10000000000 "
    cmd += " -s 0 "
    cmd += " -o " + output_filename + ".temp "
    cmd += " > " + output_filename + ".output "
    cmd += " 2> " + output_filename + ".error "
    cmd += bam_filename
    if os.path.exists( bam_filename + ".bai" ) == False:
        cmd_samtools = "samtools index " + bam_filename
        os.system( cmd_samtools )
        cmd_samtools = "samtools index -c " + bam_filename
        os.system( cmd_samtools )
    if os.path.exists( output_filename + ".temp" ) == False:
        os.system( cmd )

    fhr = open( output_filename + ".temp", "r" )
    fhw = open( output_filename, "w" )
    for line in fhr:
        chromosome, chrom_start, chrom_end, name, read_support, strand, thick_start, thick_end, item_rgb, block_count, block_sizes, block_starts = line.strip().split( "\t" )
        block_sizes_int = list( map( int, block_sizes.split( "," ) ) )
        chrom_start, chrom_end = str( int( chrom_start ) + block_sizes_int[0] + 1 ), str( int( chrom_end ) - block_sizes_int[1] )
        fhw.write( "\t".join( [chromosome, chrom_start, chrom_end, name, read_support, strand, thick_start, thick_end, item_rgb, block_count, block_sizes, block_starts] ) + "\n" )
    fhw.close()
    fhr.close()

    cmd = "rm " + output_filename + ".temp"
    os.system( cmd )


def splitTranscriptsWithQuestionableSpliceJunctions( options ):
    """
    """
    pool = multiprocessing.Pool( processes = int( options.cpu ) )
    inputfilename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_merged_transcripts.gtf"
    outputfilename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_split_transcripts_with_bad_SJ.gtf"
    if os.path.exists( outputfilename ) == True:return
    transcript_info, useless1, useless2 = readAllTranscriptsFromGTFFileInParallel( [inputfilename, "dummy", "dummy"] )
    gene_to_transcripts = {}
    for chromosome in transcript_info:
        gene_to_transcripts[chromosome] = {}
        for transcript_id in transcript_info[chromosome]:
            gene_id = ".".join( transcript_id.split( "." )[:2] )
            if gene_id not in gene_to_transcripts[chromosome]:
                gene_to_transcripts[chromosome][gene_id] = []
            gene_to_transcripts[chromosome][gene_id].append( transcript_id )

    # Generate Splice junction information using regtools
    allinputs = []
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            bam_filename = options.output_star + "/" + Run + "_for_psiclass.bam "
            output_filename = options.output_star + "/" + Run + "_SJ_regtools.bed"
            if os.path.exists( output_filename ) == False:
                allinputs.append( [bam_filename, output_filename, options] )
    pool.map( extractSJ, allinputs )

    all_SJ_info = {}
    allinputs = []
    for condition in options.mrna_md:
        all_SJ_info[condition] = {}
        for Run in options.mrna_md[condition]:
            all_SJ_info[condition][Run] = {}
            SJ_filename_regtools = options.output_star + "/" + Run + "_SJ_regtools.bed"
            allinputs.append( [SJ_filename_regtools, Run, condition, options] )
    results = pool.map( readFromRegtoolsOutput, allinputs )

    for result in results:
        SJ_info, Run, condition = result
        all_SJ_info[condition][Run] = SJ_info

    intron_info = []
    for chromosome in transcript_info:
        for transcript_id in transcript_info[chromosome]:
            # if transcript_id!="1.3.4":continue
            gene_id = ".".join( transcript_id.split( "." )[:2] )
            if len( gene_to_transcripts[chromosome][gene_id] ) == 1:continue
            if len( transcript_info[chromosome][transcript_id]["introns"] ) < 2:continue
            direction = transcript_info[chromosome][transcript_id]["direction"]
            present_in_num_conditions_direction_match = []
            present_in_num_conditions_direction_mismatch = []
            intron_coverages = []
            for condition in options.mrna_md:
                for Run in options.mrna_md[condition]:
                    intron_coverages_per_run = []
                    intron_sizes = []
                    for intron in transcript_info[chromosome][transcript_id]["introns"]:
                        intron_sizes.append( intron[1] - intron[0] + 1 )
                        if chromosome not in all_SJ_info[condition][Run]:continue
                        if str( intron[0] ) + "-" + str( intron[1] ) not in all_SJ_info[condition][Run][chromosome]:
                            intron_coverages_per_run.append( 0 )
                        else:
                            intron_coverages_per_run.append( all_SJ_info[condition][Run][chromosome][str( intron[0] ) + "-" + str( intron[1] )] )
                    if len( intron_coverages_per_run ) != 0:
                        intron_coverages.append( np.array( intron_coverages_per_run ) )
            """pprint.pprint(intron_coverages)
            sys.stdout.flush()"""
            if len( intron_coverages ) < 1:continue
            mean = list( np.mean( np.array( intron_coverages ), axis = 0 ) )
            std = []
            for ele in list( np.std( np.array( intron_coverages ), axis = 0 ) ):
                std.append( round( ele, 2 ) )
            if mean.index( min( mean ) ) != 0 and mean.index( min( mean ) ) != len( mean ) - 1:
                intron_info.append( [transcript_id, np.std( mean ), min( mean ), mean, intron_sizes] )
    all_stds = [row[1] for row in intron_info]
    all_min_mean = [row[2] for row in intron_info]
    ratios = [all_stds[i] / all_min_mean[i] for i in range( len( all_stds ) )]
    ratios = sorted( ratios )
    if len( ratios ) == 0:
        cmd = "cp " + inputfilename + " " + outputfilename
        os.system( cmd )
        return
    threshold = ratios[int( 0.75 * len( ratios ) )]
    mikado_results = {}
    """fhr=open("/work/LAS/rpwise-lab/sagnik/arabidopsis_gene_annotation/data/4_tissues_SE_PE/assemblies_psiclass_modified/combined/combined_merged_transcripts_mikado.tmap","r")
    for line in fhr:
        if line.strip().split("\t")[3] not in mikado_results:
            mikado_results[line.strip().split("\t")[3]]=[]
        mikado_results[line.strip().split("\t")[3]].append(line.strip().split("\t"))
    fhr.close()"""

    new_transcript_info = deepcopy( transcript_info )
    for row in intron_info:
        # print(row[1])
        transcript_id, std_of_intron_cov, min_mean, mean, intron_sizes = row
        # if std_of_intron_cov>=threshold and min_mean<10:
        ratio = std_of_intron_cov / min_mean
        if ratio > threshold:
            intron_indices_to_be_eliminated = [i for i in range( len( mean ) ) if mean[i] == min_mean]
            # intron_indices_to_be_eliminated.insert(0,0)
            # if len(intron_indices_to_be_eliminated)<=1:continue
            # intron_indices_to_be_eliminated.append(len(new_transcript_info[chromosome][transcript_id]["exons"])-1)

            """perfect_match="NO"
            for entry in mikado_results[transcript_id]:
                if "=" in entry[2]:
                    perfect_match="YES"
            print(" ".join([transcript_id,str(std_of_intron_cov),str(min_mean),perfect_match]),mean,intron_sizes,("large" if max(intron_sizes)>1000 else "small"))
            #print("Introns",intron_indices_to_be_eliminated,len(intron_indices_to_be_eliminated))
            """
            chromosome = transcript_id.split( "." )[0]
            previous_intron_list = new_transcript_info[chromosome][transcript_id]["introns"]
            new_exon_lists = []
            current_exon_chain = []

            for exon_num, exon in enumerate( new_transcript_info[chromosome][transcript_id]["exons"] ):
                if exon_num - 1 in intron_indices_to_be_eliminated:
                    new_exon_lists.append( current_exon_chain )
                    current_exon_chain = []
                    current_exon_chain.append( new_transcript_info[chromosome][transcript_id]["exons"][exon_num] )
                else:
                    current_exon_chain.append( new_transcript_info[chromosome][transcript_id]["exons"][exon_num] )
            new_exon_lists.append( current_exon_chain )
            """pprint.pprint(new_transcript_info[chromosome][transcript_id]["exons"])
            pprint.pprint(new_exon_lists)
            print("="*150)"""
            for exon_num, new_exon_list in enumerate( new_exon_lists ):
                new_transcript_id = transcript_id + "_isplit_" + str( exon_num )
                new_transcript_info[chromosome][new_transcript_id] = {"exons":new_exon_list,
                                                      "introns":[],
                                                      "cov":new_transcript_info[chromosome][transcript_id]["cov"],
                                                      "TPM":new_transcript_info[chromosome][transcript_id]["TPM"],
                                                      "FPKM":new_transcript_info[chromosome][transcript_id]["FPKM"],
                                                      "direction":new_transcript_info[chromosome][transcript_id]["direction"]
                                                      }
            del new_transcript_info[chromosome][transcript_id]

    writeTranscriptsToFile( [new_transcript_info, outputfilename] )
    removeRedundantTranscripts( outputfilename, outputfilename + "temp", options )
    cmd = "mv " + outputfilename + "temp " + outputfilename
    os.system( cmd )

    """cmd="mikado compare -r /work/LAS/rpwise-lab/sagnik/data/arath/transcriptome/Arabidopsis_thaliana.TAIR10.43.modified.gtf "
    cmd+=" -p "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_split_transcripts_with_bad_SJ.gtf"
    cmd+=" -o "+options.output_assemblies_psiclass_terminal_exon_length_modified+"/combined/combined_split_transcripts_with_bad_SJ_mikado "
    os.system(cmd)"""

