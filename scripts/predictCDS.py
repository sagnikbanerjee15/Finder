
from scripts.fileReadWriteOperations import *
import os


def convertToGenomicCoordinate( transcriptomic_coordinate, exon_list_genomic, transcript_id ):
    """
    """
    exon_list_transcriptomic = []
    for exon in exon_list_genomic:
        exon_start, exon_end = exon
        exon_length = exon_end - exon_start + 1
        if len( exon_list_transcriptomic ) == 0:
            exon_list_transcriptomic.append( [1, exon_length] )
        else:
            exon_list_transcriptomic.append( [exon_list_transcriptomic[-1][1] + 1, exon_length + exon_list_transcriptomic[-1][1]] )

    """if transcript_id=="1.62.0":
        print(transcriptomic_coordinate)"""
    for exon_num, exon in enumerate( exon_list_transcriptomic ):
        exon_start, exon_end = exon
        if exon_start <= transcriptomic_coordinate <= exon_end:
            """if transcript_id=="1.62.0":
                print(exon_list_genomic[exon_num][0],transcriptomic_coordinate,exon_start,exon_list_genomic[exon_num][0]+transcriptomic_coordinate-exon_start)"""
            return exon_list_genomic[exon_num][0] + transcriptomic_coordinate - exon_start

    """print(exon_list_transcriptomic,transcriptomic_coordinate)
    print("="*100)"""


def findCDS( options ):
    """
    """
    gtf_filename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_split_transcripts_with_bad_SJ_redundancy_removed.gtf"
    combined_ultra_long_introns_redundancy_removed = readAllTranscriptsFromGTFFileInParallel( [gtf_filename, "dummy", "dummy"] )[0]
    output_gtf_filename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.gtf"
    if os.path.exists( output_gtf_filename ) == True:return
    for transcript_id in combined_ultra_long_introns_redundancy_removed:
        if len( combined_ultra_long_introns_redundancy_removed[transcript_id]["exons"] ) == 1:
            combined_ultra_long_introns_redundancy_removed[transcript_id]["direction"] = "."

    writeTranscriptsToFile( [combined_ultra_long_introns_redundancy_removed, gtf_filename, 0] )

    os.chdir( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined" )
    # Convert from gtf to fasta
    cmd = "gffread "
    cmd += " -w " + gtf_filename[:-3] + "fasta "
    cmd += " -g " + options.genome
    cmd += " " + gtf_filename
    os.system( cmd )

    cmd = f" codan.py "
    cmd += f" -t {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/combined_ultra_long_introns_redundancy_removed_genemarkST_output.gff3 "
    cmd += f" -m /softwares/CODAN/CodAn-1.2/models/${options.organism_model}_full "
    cmd += f" -c {options.cpu} "
    cmd += f" -o {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/cds_predict "
    cmd += f" 1> {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/cds_predict.output "
    cmd += f" 2> {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/cds_predict.error "
    if os.path.exists( f"{options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/cds_predict/annotation.gtf" ) == False:
        os.system( cmd )

    transcript_info = readAllTranscriptsFromGTFFileInParallel( [gtf_filename, "dummy", "dummy"] )[0]
    transcript_to_CDS_start_end = {}
    fhr = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/cds_predict/annotation.gtf", "r" )
    for line in fhr:
        if line.strip().split( "\t" )[2] != "CDS":continue
        # print(line.strip().split("\t"))
        transcript_id = line.strip().split( "\t" )[0]
        start, end = int( line.strip().split( "\t" )[3] ), int( line.strip().split( "\t" )[4] )
        transcript_to_CDS_start_end[transcript_id] = [start, end]
    fhr.close()

    output_filename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.gtf"
    fhw = open( output_filename, "w" )

    for transcript_id in transcript_info:
        chromosome = transcript_info[transcript_id]["chromosome"]
        transcript_start = transcript_info[transcript_id]["exons"][0][0]
        transcript_end = transcript_info[transcript_id]["exons"][-1][1]
        annotator = transcript_info[transcript_id]["annotator"]
        gene_id = ".".join( transcript_id.split( "." )[:2] )
        desc = "gene_id \"" + gene_id + "\"; transcript_id \"" + transcript_id + "\"; "
        desc += "FPKM \"" + str( transcript_info[transcript_id]["FPKM"] ) + "\"; "
        desc += "TPM \"" + str( transcript_info[transcript_id]["TPM"] ) + "\"; "
        desc += "cov \"" + str( transcript_info[transcript_id]["cov"] ) + "\"; "
        if transcript_id not in transcript_to_CDS_start_end:
            row = [chromosome,
             annotator,
             "transcript",
             str( transcript_start ),
             str( transcript_end ),
             "1000",
             transcript_info[transcript_id]["direction"],
             ".",
             desc
             ]
            fhw.write( "\t".join( row ) + "\n" )

            for exon in transcript_info[transcript_id]["exons"]:
                exon_start, exon_end = exon[0], exon[1]
                row = [chromosome,
                 annotator,
                 "exon",
                 str( exon_start ),
                 str( exon_end ),
                 "1000",
                 transcript_info[transcript_id]["direction"],
                 ".",
                 desc
                 ]
                fhw.write( "\t".join( row ) + "\n" )
        else:
            CDS_start, CDS_end = transcript_to_CDS_start_end[transcript_id]
            # CDS_start_genome_reference=transcript_start+CDS_start-1
            # CDS_end_genome_reference=transcript_start+CDS_end-1
            if transcript_info[transcript_id]["direction"] == "+" or len( transcript_info[transcript_id]["exons"] ) == 1 or transcript_info[transcript_id]["direction"] == ".":
                CDS_start_genome_reference = convertToGenomicCoordinate( CDS_start, transcript_info[transcript_id]["exons"], transcript_id )
                CDS_end_genome_reference = convertToGenomicCoordinate( CDS_end, transcript_info[transcript_id]["exons"], transcript_id )
            elif transcript_info[transcript_id]["direction"] == "-":
                """if "3.15366.0"==transcript_id:
                    print("doing")
                    print(CDS_start,CDS_end)"""
                transcript_length = sum( [exon[1] - exon[0] + 1 for exon in transcript_info[transcript_id]["exons"]] )
                # print(transcript_length,CDS_start,CDS_end,transcript_length-CDS_end+1,transcript_length-CDS_start+1)
                CDS_start_genome_reference = convertToGenomicCoordinate( transcript_length - CDS_end + 1, transcript_info[transcript_id]["exons"], transcript_id )
                CDS_end_genome_reference = convertToGenomicCoordinate( transcript_length - CDS_start + 1, transcript_info[transcript_id]["exons"], transcript_id )
                """if "3.15366.0"==transcript_id:
                    print(CDS_start_genome_reference,CDS_end_genome_reference)"""

            row = [chromosome,
             annotator,
             "transcript",
             str( transcript_start ),
             str( transcript_end ),
             "1000",
             transcript_info[transcript_id]["direction"],
             ".",
             desc
             ]
            fhw.write( "\t".join( row ) + "\n" )

            for exon in transcript_info[transcript_id]["exons"]:
                exon_start, exon_end = exon[0], exon[1]
                row = [chromosome,
                 annotator,
                 "exon",
                 str( exon_start ),
                 str( exon_end ),
                 "1000",
                 transcript_info[transcript_id]["direction"],
                 ".",
                 desc
                 ]
                fhw.write( "\t".join( row ) + "\n" )

            if transcript_info[transcript_id]["direction"] != ".":
                start_writing = 0
                for exon in transcript_info[transcript_id]["exons"]:
                    exon_start, exon_end = exon[0], exon[1]
                    print( transcript_id, exon_start, CDS_start_genome_reference, exon_end, CDS_end_genome_reference )
                    if start_writing == 0 and exon_start <= CDS_start_genome_reference <= exon_end and ( exon_start <= CDS_end_genome_reference <= exon_end ) == False:
                        start_writing = 1
                        row = [chromosome,
                         annotator,
                         "CDS",
                         str( CDS_start_genome_reference ),
                         str( exon_end ),
                         "1000",
                         transcript_info[transcript_id]["direction"],
                         ".",
                         desc
                         ]
                        fhw.write( "\t".join( row ) + "\n" )
                    elif start_writing == 0 and exon_start <= CDS_start_genome_reference <= exon_end and exon_start <= CDS_end_genome_reference <= exon_end:
                        start_writing = 1
                        row = [chromosome,
                         annotator,
                         "CDS",
                         str( CDS_start_genome_reference ),
                         str( CDS_end_genome_reference ),
                         "1000",
                         transcript_info[transcript_id]["direction"],
                         ".",
                         desc
                         ]
                        fhw.write( "\t".join( row ) + "\n" )

                    elif start_writing == 1 and CDS_start_genome_reference <= exon_start and CDS_end_genome_reference >= exon_end:
                        row = [chromosome,
                         annotator,
                         "CDS",
                         str( exon_start ),
                         str( exon_end ),
                         "1000",
                         transcript_info[transcript_id]["direction"],
                         ".",
                         desc
                         ]
                        fhw.write( "\t".join( row ) + "\n" )
                    elif start_writing == 1 and exon_start <= CDS_end_genome_reference <= exon_end:
                        start_writing = 0
                        row = [chromosome,
                         annotator,
                         "CDS",
                         str( exon_start ),
                         str( CDS_end_genome_reference ),
                         "1000",
                         transcript_info[transcript_id]["direction"],
                         ".",
                         desc
                         ]
                        fhw.write( "\t".join( row ) + "\n" )
            else:
                pass

    fhw.close()

    # Adding CDS frame info by calling readAllTranscriptsFromGTFFileInParallel
    all_info_with_CDS = readAllTranscriptsFromGTFFileInParallel( [options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.gtf", "dummy", "dummy"] )[0]
    fhw = open( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_temp.gtf", "w" )

    for transcript_id in all_info_with_CDS:
        chromosome = all_info_with_CDS[transcript_id]["chromosome"]
        this_transcript = all_info_with_CDS[transcript_id]
        transcript_start = this_transcript["exons"][0][0]
        transcript_end = this_transcript["exons"][-1][1]
        gene_id = ".".join( transcript_id.split( "." )[:2] )
        desc = "gene_id \"" + gene_id + "\"; transcript_id \"" + transcript_id + "\"; "
        desc += "FPKM \"" + str( this_transcript["FPKM"] ) + "\"; "
        desc += "TPM \"" + str( this_transcript["TPM"] ) + "\"; "
        desc += "cov \"" + str( this_transcript["cov"] ) + "\"; "

        row = [chromosome,
             annotator,
             "transcript",
             str( transcript_start ),
             str( transcript_end ),
             "1000",
             this_transcript["direction"],
             ".",
             desc
             ]
        fhw.write( "\t".join( row ) + "\n" )

        for exon in this_transcript["exons"]:
            exon_start, exon_end = exon[0], exon[1]
            row = [chromosome,
             annotator,
             "exon",
             str( exon_start ),
             str( exon_end ),
             "1000",
             this_transcript["direction"],
             ".",
             desc
             ]
            fhw.write( "\t".join( row ) + "\n" )

        for cds_num, cds in enumerate( this_transcript["cds"] ):
            cds_start, cds_end = cds[0], cds[1]
            row = [chromosome,
             annotator,
             "CDS",
             str( cds_start ),
             str( cds_end ),
             "1000",
             this_transcript["direction"],
             str( this_transcript["cds_frame"][cds_num] ),
             desc
             ]
            fhw.write( "\t".join( row ) + "\n" )
    fhw.close()
    cmd = "mv " + options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_temp.gtf "
    cmd += options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.gtf"
    os.system( cmd )

