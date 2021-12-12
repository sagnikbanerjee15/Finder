#! /usr/bin/env python3
"""
Load these modules

module purge
module load bedtools2
module load samtools
module load cufflinks
module load bedops
module load genometools
module load python
module load py-pandas/0.23.4-py3-sx6iffy
"""

import argparse
import os
import parser
import pickle
import subprocess
import sys

import numpy as np
import pandas as pd


def extractExons( gtffilename, outputprefix ):
    """
    """
    cmd = "python hisat2_extract_exons.py " + gtffilename + " > " + outputprefix + "_exons"
    os.system( cmd )


def readExons( exonfilename ):
    """
    """
    exons = {}
    fhr = open( exonfilename, "r" )
    for line in fhr:
        chromosome, start, end, direction = line.strip().split( "\t" )
        start, end = int( start ) + 1, int( end ) + 1
        if chromosome not in exons:
            exons[chromosome] = []
        exons[chromosome].append( [start, end] )
    fhr.close()
    return exons


def readGenomicCounts( bedfilename ):
    """
    """
    coverage_info = {}
    fhr = open( bedfilename, "r" )
    for line_num, line in enumerate( fhr ):
        chromosome, start, end, cov = line.strip().split( "\t" )
        start, end, cov = int( start ) + 1, int( end ), int( float( cov ) )
        if cov == 0:continue
        if chromosome not in coverage_info:
            coverage_info[chromosome] = {}
        for i in range( start, end + 1 ):
            coverage_info[chromosome][i] = cov
        if line_num % 10000 == 0:
            print( line_num, bedfilename )
            sys.stdout.flush()
    fhr.close()
    return coverage_info


def transferCounts( coverage_info, transcript_info, transcriptome_out_bed ):
    """
    """
    for transcript_id in transcript_info:
        for eachexon in transcript_info[transcript_id]["exons"]:
            start, end = int( eachexon.split( "_" )[0] ), int( eachexon.split( "_" )[1] )
            exon_coverage = []
            for i in range( start, end + 1 ):
                try:
                    exon_coverage.append( coverage_info[transcript_info[transcript_id]["chromosome"]][i] )
                except KeyError:
                    # print(transcript_id,"KEYERROR",i,start,end)
                    exon_coverage.append( 0 )
            transcript_info[transcript_id]["bed_cov"].append( exon_coverage )
    return transcript_info


def readAllTranscriptsFromGTFFileInParallel( gtf_filename ):
    whole_annotations = {}
    fhr = open( gtf_filename, "r" )
    for line in fhr:
        if line[0] == "#":continue
        chromosome, psiclass, structure, start, end, useless1, direction, useless2, desc = line.strip().split( "\t" )
        exon = start + "-" + end
        tpm = 0
        fpkm = 0
        cov = 0
        for ele in desc.split( ";" ):
            if "gene_id" in ele:
                gene_id = ele.split()[-1].strip( "\"" )
            if "transcript_id" in ele:
                transcript_id = ele.split()[-1].strip( "\"" )
            try:
                if "cov" in ele:
                    cov = float( ele.split()[-1].strip( "\"" ) )
                if "TPM" in ele:
                    tpm = float( ele.split()[-1].strip( "\"" ) )
                if "FPKM" in ele:
                    fpkm = float( ele.split()[-1].strip( "\"" ) )
            except ValueError:
                tpm = 0
                fpkm = 0
                cov = 0
        if structure == "transcript":
            continue
        if chromosome not in whole_annotations:
            whole_annotations[chromosome] = {}

        if transcript_id not in whole_annotations[chromosome]:
            whole_annotations[chromosome][transcript_id] = {"exons":[],
                                                        "introns":[],
                                                        "cds":[],
                                                        "cds_frame":[],
                                                        "direction":direction,
                                                        "TPM":tpm,
                                                        "cov":cov,
                                                        "gene_id":gene_id,
                                                        "FPKM":fpkm,
                                                        "transcript_start":0,
                                                        "transcript_end":0
                                                        }
        if structure == "exon":
            whole_annotations[chromosome][transcript_id]["exons"].append( [int( start ), int( end )] )
        elif structure == "CDS":
            whole_annotations[chromosome][transcript_id]["cds"].append( [int( start ), int( end )] )

        # Construct the introns for each transcript
        if len( whole_annotations[chromosome][transcript_id]["exons"] ) > 1:
            exon1 = int( whole_annotations[chromosome][transcript_id]["exons"][-2][1] )
            exon2 = int( whole_annotations[chromosome][transcript_id]["exons"][-1][0] )
            whole_annotations[chromosome][transcript_id]["introns"].append( [exon1 + 1, exon2 - 1] )
    fhr.close()

    for chromosome in whole_annotations:
        for transcript_id in whole_annotations[chromosome]:
            whole_annotations[chromosome][transcript_id]["transcript_start"] = whole_annotations[chromosome][transcript_id]["exons"][0][0]
            whole_annotations[chromosome][transcript_id]["transcript_end"] = whole_annotations[chromosome][transcript_id]["exons"][-1][1]
        for transcript_id in whole_annotations[chromosome]:
            if whole_annotations[chromosome][transcript_id]["direction"] == "+":
                whole_annotations[chromosome][transcript_id]["cds_frame"].append( 0 )
                for cds_num, cds in enumerate( whole_annotations[chromosome][transcript_id]["cds"] ):
                    if cds_num == 0:continue
                    length_of_CDS_sequence = sum( [cds[1] - cds[0] + 1 for cds in whole_annotations[chromosome][transcript_id]["cds"][:cds_num]] )
                    if length_of_CDS_sequence % 3 == 0:
                        whole_annotations[chromosome][transcript_id]["cds_frame"].append( 0 )
                    elif length_of_CDS_sequence % 3 == 1:
                        whole_annotations[chromosome][transcript_id]["cds_frame"].append( 2 )
                    elif length_of_CDS_sequence % 3 == 2:
                        whole_annotations[chromosome][transcript_id]["cds_frame"].append( 1 )
            elif whole_annotations[chromosome][transcript_id]["direction"] == "-":
                whole_annotations[chromosome][transcript_id]["cds_frame"].append( 0 )
                for cds_num, cds in enumerate( whole_annotations[chromosome][transcript_id]["cds"][::-1] ):
                    if cds_num == 0:continue
                    length_of_CDS_sequence = sum( [cds[1] - cds[0] + 1 for cds in whole_annotations[chromosome][transcript_id]["cds"][::-1][:cds_num]] )
                    if length_of_CDS_sequence % 3 == 0:
                        whole_annotations[chromosome][transcript_id]["cds_frame"].append( 0 )
                    elif length_of_CDS_sequence % 3 == 1:
                        whole_annotations[chromosome][transcript_id]["cds_frame"].append( 2 )
                    elif length_of_CDS_sequence % 3 == 2:
                        whole_annotations[chromosome][transcript_id]["cds_frame"].append( 1 )
                rev = whole_annotations[chromosome][transcript_id]["cds_frame"][::-1]
                whole_annotations[chromosome][transcript_id]["cds_frame"] = rev
            else:
                pass
    return whole_annotations


def writeTranscriptsToFile( transcript_info, outputfilename ):
    """
    """
    fhw = open( outputfilename, "w" )
    for chromosome in transcript_info:
        for transcript_id in transcript_info[chromosome]:
            transcript_start = transcript_info[chromosome][transcript_id]["exons"][0][0]
            transcript_end = transcript_info[chromosome][transcript_id]["exons"][-1][1]
            gene_id = ".".join( transcript_id.split( "." )[:2] )
            desc = "gene_id \"" + gene_id + "\"; transcript_id \"" + transcript_id + "\"; "
            desc += "FPKM \"" + str( transcript_info[chromosome][transcript_id]["FPKM"] ) + "\"; "
            desc += "TPM \"" + str( transcript_info[chromosome][transcript_id]["TPM"] ) + "\"; "
            desc += "cov \"" + str( transcript_info[chromosome][transcript_id]["cov"] ) + "\"; "
            row = [chromosome,
             "PsiCLASS",
             "transcript",
             str( transcript_start ),
             str( transcript_end ),
             "1000",
             transcript_info[chromosome][transcript_id]["direction"],
             ".",
             desc
             ]
            fhw.write( "\t".join( row ) + "\n" )
            for exon in transcript_info[chromosome][transcript_id]["exons"]:
                exon_start, exon_end = exon[0], exon[1]
                row = [chromosome,
                 "PsiCLASS",
                 "exon",
                 str( exon_start ),
                 str( exon_end ),
                 "1000",
                 transcript_info[chromosome][transcript_id]["direction"],
                 ".",
                 desc
                 ]
                fhw.write( "\t".join( row ) + "\n" )

            for intron in transcript_info[chromosome][transcript_id]["introns"]:
                intron_start, intron_end = intron[0], intron[1]
                row = [chromosome,
                 "PsiCLASS",
                 "intron",
                 str( intron_start ),
                 str( intron_end ),
                 "1000",
                 transcript_info[chromosome][transcript_id]["direction"],
                 ".",
                 desc
                 ]
                fhw.write( "\t".join( row ) + "\n" )

            if "cds" in transcript_info[chromosome][transcript_id]:
                for cds in transcript_info[chromosome][transcript_id]["cds"]:
                    exon_start, exon_end = cds[0], cds[1]
                    row = [chromosome,
                     "PsiCLASS",
                     "CDS",
                     str( exon_start ),
                     str( exon_end ),
                     "1000",
                     transcript_info[chromosome][transcript_id]["direction"],
                     ".",
                     desc
                     ]
                    fhw.write( "\t".join( row ) + "\n" )
    fhw.close()


def readTranscriptInfoFromGTFFile( gtffilename, exons_overlapping_with_introns_bedfilename, portion_of_exons_overlapping_with_introns_bedfilename, options ):
    """
    """

    """print("Entering readTranscriptInfoFromGTFFile")
    sys.stdout.flush()"""
    if os.path.exists( exons_overlapping_with_introns_bedfilename ) == False:
        cmd = "bedops " + " --element-of 1 "
        cmd += " <(" + "gff2bed " + " < <(" + "gffread " + " -E " + gtffilename + " -o-|tail -n +2) |grep exon) "
        cmd += " <(" + "bedops " + " --intersect <(" + "gff2bed " + " < <(gt gff3 -addids yes -addintrons yes "
        cmd += " <(" + "gffread " + " -E " + gtffilename + " -o-|tail -n +3) |grep exon)) "
        cmd += " <(" + "gff2bed " + " < <(gt gff3 -addids yes -addintrons yes "
        cmd += " <(" + "gffread " + " -E " + gtffilename + " -o-|tail -n +3) |grep intron))) "
        cmd += " > " + exons_overlapping_with_introns_bedfilename
        # subprocess.check_call(['bash', '-c', cmd])
        # print(cmd)
        writeTranscriptsToFile( readAllTranscriptsFromGTFFileInParallel( gtffilename ), gtffilename[:-4] + "_appended_info.gtf" )

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
        # print(cmd)
        os.system( cmd )
        sys.stdout.flush()
        # sys.exit()
        # print(cmd)

    overlaps = {}
    exons_to_transcripts = {}
    fhr = open( exons_overlapping_with_introns_bedfilename, "r" )
    for line in fhr:
        chromosome, start, end = line.strip().split( "\t" )[:3]
        if chromosome not in overlaps:
            overlaps[chromosome] = []
        overlaps[chromosome].append( str( int( start ) + 1 ) + "_" + str( end ) )
    fhr.close()

    """portion_of_overlaps={}
    fhr=open(portion_of_exons_overlapping_with_introns_bedfilename,"r")
    for line in fhr:
        chromosome,start,end=line.strip().split("\t")[:3]
        if chromosome not in portion_of_overlaps:
            portion_of_overlaps[chromosome]=[]
        portion_of_overlaps[chromosome].append([int(start)+1,int(end)])
    fhr.close()
    for chromosome in portion_of_overlaps:
        portion_of_overlaps[chromosome]=np.array(portion_of_overlaps[chromosome])
        portion_of_overlaps[chromosome]=portion_of_overlaps[chromosome][np.argsort(portion_of_overlaps[chromosome][:, 0])]
        #portion_of_overlaps[chromosome].columns=["start","end"]"""

    all_exons_in_introns = {}
    iterator = 0
    transcript_info = {}
    fhr = open( gtffilename, "r" )
    """print(gtffilename)
    print("NODE_1_length_3753_cov_430.083741_g0_i0.mrna1_1" in open(gtffilename,"r").read())"""
    for line in fhr:
        chromosome, method, structure, start, end, useless1, direction, useless2, info = line.strip().split( "\t" )
        if chromosome not in all_exons_in_introns:
            all_exons_in_introns[chromosome] = {}
            iterator = 0
        if structure == "transcript":
            for ele in info.split( ";" ):
                if "gene_id" in ele:
                    gene = ele.split()[-1].strip( "\"" )
                if "transcript_id" in ele:
                    transcript = ele.split()[-1].strip( "\"" )
            gene = gene.split()[-1].strip( "\"" )
            transcript = transcript.split()[-1].strip( "\"" )
            RPKM = -1
            cov = -1
            transcript_info[transcript] = {"gene":gene, "RPKM":RPKM, "cov_gtf":cov, "exons":[], "bed_cov":[], "chromosome":chromosome, "exons_in_introns":[], "introns":[], "direction":direction}
        elif structure == "exon":
            for ele in info.split( ";" ):
                if "gene_id" in ele:
                    gene_id = ele.split()[-1].strip( "\"" )
                if "transcript_id" in ele:
                    transcript = ele.split()[-1].strip( "\"" )
            RPKM = -1
            cov = -1
            """if transcript=="NODE_1_length_3753_cov_430.083741_g0_i0.mrna1_1":
                print("Present")
                sys.stdout.flush()"""
            if transcript not in transcript_info:
                transcript_info[transcript] = {"gene":gene, "RPKM":RPKM, "cov_gtf":cov, "exons":[], "bed_cov":[], "chromosome":chromosome, "exons_in_introns":[], "introns":[], "direction":direction}
            transcript_info[transcript]["exons"].append( str( start ) + "_" + str( end ) )
            if len( transcript_info[transcript]["exons"] ) > 1:
                exon1 = int( transcript_info[transcript]["exons"][-2].split( "_" )[1] )
                exon2 = int( transcript_info[transcript]["exons"][-1].split( "_" )[0] )
                transcript_info[transcript]["introns"].append( [exon1 + 1, exon2 - 1] )
            if transcript_info[transcript]["chromosome"] in overlaps and str( start ) + "_" + str( end ) in overlaps[transcript_info[transcript]["chromosome"]]:
                transcript_info[transcript]["exons_in_introns"].append( str( start ) + "_" + str( end ) )
                """print(chromosome+":"+str(start)+"-"+str(end))
                if str(start)+"_"+str(end) not in all_exons_in_introns[chromosome]:
                    x=iterator
                    while x<len(portion_of_overlaps[chromosome]):
                        overlap_start,overlap_end=portion_of_overlaps[chromosome][x][0],portion_of_overlaps[chromosome][x][1]
                        #print(chromosome+":"+str(start)+"-"+str(end),~((overlap_start>=int(end) and overlap_end>=int(end)) or (overlap_start<=int(start) and overlap_end<=int(start))))
                        if ~((overlap_start>=int(end) and overlap_end>=int(end)) or (overlap_start<=int(start) and overlap_end<=int(start))):
                            print(chromosome+":"+str(start)+"-"+str(end),chromosome+":"+str(overlap_start)+"-"+str(overlap_end))
                        else:
                            pass
                        if int(end)<overlap_start:
                            break
                        x+=1
                    iterator=x-10 if x-10>=0 else 0
                    all_exons_in_introns[chromosome][str(start)+"_"+str(end)]=x"""
                # print(chromosome, portion_of_overlaps[chromosome].loc[~(portion_of_overlaps[chromosome]['start'] > int(end)) | (portion_of_overlaps[chromosome]['end'] < int(start))].to_string())
                """print(chromosome, portion_of_overlaps[chromosome].loc[((portion_of_overlaps[chromosome]['start'] > int(start)) & (portion_of_overlaps[chromosome]['end'] < int(end)) |
                                                                       (portion_of_overlaps[chromosome]['start'] > int(end)) & (portion_of_overlaps[chromosome]['end'] < int(start)) |
                                                                       (portion_of_overlaps[chromosome]['start'] > int(end)) & (portion_of_overlaps[chromosome]['end'] < int(start)) |
                                                                       (portion_of_overlaps[chromosome]['start'] > int(end)) & (portion_of_overlaps[chromosome]['end'] < int(start))
                                                                       )].to_string())
                """
            if chromosome + "_" + str( start ) + "_" + str( end ) not in exons_to_transcripts:
                exons_to_transcripts[chromosome + "_" + str( start ) + "_" + str( end )] = []
            exons_to_transcripts[chromosome + "_" + str( start ) + "_" + str( end )].append( transcript )
    fhr.close()
    return [transcript_info, exons_to_transcripts]


def transferCountsFromGenomeToExons( bedfilename, exonfilename, gtffilename, transcriptome_out_bed, exons_overlapping_with_introns_bedfilename, portion_of_exons_overlapping_with_introns_bedfilename, pklfile, options ):
    """
    """
    transcript, exons_to_transcripts = readTranscriptInfoFromGTFFile( gtffilename, exons_overlapping_with_introns_bedfilename, portion_of_exons_overlapping_with_introns_bedfilename, options )
    # exons=readExons(exonfilename)
    coverage_info = readGenomicCounts( bedfilename )
    # print("Coverage info read ",bedfilename)
    sys.stdout.flush()
    transcript = transferCounts( coverage_info, transcript, transcriptome_out_bed )
    # print("Gene counts transferred",bedfilename)
    sys.stdout.flush()
    pickle.dump( [transcript, exons_to_transcripts], open( pklfile, "wb" ) )
    """for transcript_id in transcript:
        #print(transcript_id,transcript[transcript_id])
        for num,exon in enumerate(transcript[transcript_id]["exons"]):
            print(("*" if exon in transcript[transcript_id]["exons_in_introns"] else ""),transcript_id,exon,transcript[transcript_id]["bed_cov"][num])"""
    # os.system("rm "+exons_overlapping_with_introns_bedfilename)


def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python transfer_genomic_nucl_counts_to_transcriptome.py --help for more details
    """
    parser = argparse.ArgumentParser( prog = "transfer_genomic_nucl_counts_to_transcriptome.py", description = "Uses aligned RNA-Seq data to derive read counts on genomic scale per nucleotide" )
    mutex_parser = parser.add_mutually_exclusive_group()
    mutex_parser.add_argument( "--bamfile", "-b", help = "Enter the name of the bamfile containing all the alignments. The program will sort the file if unsorted file is provided." )
    mutex_parser.add_argument( "--bedfile", "-bd", help = "Enter the name of the bedfile containing the per nucleotide coverage on the genome locus. You can generate this file by issuing the command: bedtools genomecov -bga -split -ibam <bamfilename> > <bedfile>" )
    parser.add_argument( "--cpu", "-n", help = "Enter the number of CPUs to be used", default = 1 )
    parser.add_argument( "--gtffilename", "-g", help = "Enter the name of the gtf file containing all the transcript sequences. Please note you do not need CDS annotation." )
    parser.add_argument( "--outprefix", "-p", help = "Enter the output prefix", required = True )
    parser.add_argument( "--force", "-f", help = "Set it enforce recalculation of all the steps", action = "store_true" )
    # Suppressed arguments
    parser.add_argument( "--softwares", help = argparse.SUPPRESS )  # Sets up the full path to the dependent softwares like STAR, Spring, etc

    return parser.parse_args()


def generatebedfile( options ):
    """
    Generates bedfile if only bamfile is provided
    """
    if options.bamfile is not None:
        # Check if the bamfile is sorted
        cmd = "samtools " + " view -H " + options.bamfile + " | grep @HD|grep coordinate > " + options.outprefix + "_check_for_sorted"
        os.system( cmd )

        # Sort the file if bamfile is not sorted
        if os.stat( options.outprefix + "_check_for_sorted" ).st_size == 0:
            cmd = "samtools " + " sort -@ " + str( options.cpu ) + " -o " + options.outprefix + "_sorted.bam " + options.bamfile
            os.system( cmd )
            options.bamfile = options.outprefix + "_sorted.bam"
        os.system( "rm " + options.outprefix + "_check_for_sorted" )

        # Generate bedfile
        options.bedfile = options.outprefix + "_genome_cov.bed"
        if os.path.exists( options.bedfile ) == False or options.force == True:
            cmd = "bedtools " + " genomecov -bga -split -ibam " + options.bamfile + " > " + options.bedfile
            os.system( cmd )
    return options


def main():
    commandLineArg = sys.argv
    if len( commandLineArg ) == 1:
        print( "Please use the --help option to get usage information" )
    options = parseCommandLineArguments()

    ########################################################################################################
    # Setting up links to the softwares
    ########################################################################################################
    options.softwares = {}
    cmd = "whereis transfer_genomic_nucl_counts_to_transcriptome.py > " + options.outprefix + "_finding_program"
    os.system( cmd )
    program_directory = open( options.outprefix + "_finding_program", "r" ).read().split( ":" )[-1].strip()
    if program_directory != "":
        # PATH NOT set to downloaded directory
        lib_path = "/".join( program_directory.split( "/" )[:-1] ) + "/lib"
    else:
        # PATH IS set to downloaded directory
        if "/" in sys.argv[0]:
            lib_path = "/".join( sys.argv[0].split( "/" )[:-1] ) + "/lib"
        else:
            lib_path = "lib"

    os.system( "rm " + options.outprefix + "_finding_program" )

    options = generatebedfile( options )
    # extractExons(options.gtffilename,options.outprefix)
    transferCountsFromGenomeToExons( options.bedfile,
                                    options.outprefix + "_exons",
                                    options.gtffilename,
                                    options.outprefix + "_transcriptome_cov.bed",
                                    options.gtffilename[:-4] + "_exons_overlapping_with_introns.bed",
                                    options.outprefix + "_portion_of_exons_overlapping_with_introns.bed",
                                    options.outprefix + "_all_info.pkl", options )


if __name__ == "__main__":
    main()
