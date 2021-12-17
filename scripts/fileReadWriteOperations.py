#! /usr/bin/env python3

###################################################################################
# File consists of functions involving read write operations
###################################################################################

import os
import time
import glob


def divide_chunks( l, n ):
    # looping till length l
    for i in range( 0, len( l ), n ):
        yield l[i:i + n]


def samtoolsQuickCheck( filename, options ):
    cmd = "samtools " + " quickcheck "
    cmd += filename + " 2> "
    cmd += filename + ".quickcheck "
    os.system( cmd )
    time.sleep( 2 )
    if os.stat( filename ).st_size == 0:
        return 0
    if os.stat( filename + ".quickcheck" ).st_size != 0:
        os.system( "rm " + filename + ".quickcheck" )
        return 1
    os.system( "rm " + filename + ".quickcheck" )
    return 0


def isValidLocation( location ):
    """
    Checks if the location provided is valid
    """
    if ( ".." in location ):
        return 1
    elif ( "/" in location ):
        return 1
    elif ( len( set( location ) ) == 1 and location[0] == ' ' ):
        return 0
    elif ( len( set( location ) ) == 0 ):
        return 0


def expandGzippedFiles( options, logger_proxy, logging_mutex ):
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            if options.mrna_md[condition][Run]["downloaded_from_NCBI"] == 0:
                location = options.mrna_md[condition][Run]["location_directory"]
                if os.path.exists( f"{location}/{Run}.fastq.gz" ) == True:
                    cmd = f"gunzip -c "
                    cmd += f"{location}/{Run}.fastq.gz "
                    cmd += " > "
                    cmd += f"{options.raw_data_downloaded_from_NCBI}/{Run}.fastq"
                    with logging_mutex:
                        logger_proxy.info( "Fastq data is expanded" )
                    os.system( cmd )
                    options.mrna_md[condition][Run]["location_directory"] = options.raw_data_downloaded_from_NCBI
                    options.mrna_md[condition][Run]["downloaded_from_NCBI"] = 1
                elif os.path.exists( f"{location}/{Run}.fq.gz" ) == True:
                    cmd = f"gunzip -c "
                    cmd += f"{location}/{Run}.fq.gz "
                    cmd += " > "
                    cmd += f"{options.raw_data_downloaded_from_NCBI}/{Run}.fastq"
                    with logging_mutex:
                        logger_proxy.info( "Fastq data is expanded" )
                    os.system( cmd )
                    options.mrna_md[condition][Run]["location_directory"] = options.raw_data_downloaded_from_NCBI
                    options.mrna_md[condition][Run]["downloaded_from_NCBI"] = 1
                elif os.path.exists( f"{location}/{Run}_1.fastq.gz" ) == True and os.path.exists( f"{location}/{Run}_2.fastq.gz" ) == True:
                    cmd = f"gunzip -c "
                    cmd += f"{location}/{Run}_1.fastq.gz "
                    cmd += " > "
                    cmd += f"{options.raw_data_downloaded_from_NCBI}/{Run}_1.fastq"
                    with logging_mutex:
                        logger_proxy.info( "Fastq data is expanded" )
                    os.system( cmd )

                    cmd = f"gunzip -c "
                    cmd += f"{location}/{Run}_2.fastq.gz "
                    cmd += " > "
                    cmd += f"{options.raw_data_downloaded_from_NCBI}/{Run}_2.fastq"
                    with logging_mutex:
                        logger_proxy.info( "Fastq data is expanded" )
                    os.system( cmd )
                    options.mrna_md[condition][Run]["location_directory"] = options.raw_data_downloaded_from_NCBI
                    options.mrna_md[condition][Run]["downloaded_from_NCBI"] = 1
                elif os.path.exists( f"{location}/{Run}_1.fq.gz" ) == True and os.path.exists( f"{location}/{Run}_2.fq.gz" ) == True:
                    cmd = f"gunzip -c "
                    cmd += f"{location}/{Run}_1.fq.gz "
                    cmd += " > "
                    cmd += f"{options.raw_data_downloaded_from_NCBI}/{Run}_1.fastq"
                    with logging_mutex:
                        logger_proxy.info( "Fastq data is expanded" )
                    os.system( cmd )

                    cmd = f"gunzip -c "
                    cmd += f"{location}/{Run}_2.fq.gz "
                    cmd += " > "
                    cmd += f"{options.raw_data_downloaded_from_NCBI}/{Run}_2.fastq"
                    with logging_mutex:
                        logger_proxy.info( "Fastq data is expanded" )
                    os.system( cmd )
                    options.mrna_md[condition][Run]["location_directory"] = options.raw_data_downloaded_from_NCBI
                    options.mrna_md[condition][Run]["downloaded_from_NCBI"] = 1


def readMetaDataFile( options, logger_proxy, logging_mutex ):
    all_samples = {}
    small_rna_samples = {}
    fhr = open( options.metadatafile, "r" )
    for line in fhr:
        if "BioProject" in line:continue
        BioProject, Run, condition, desc, Date, read_length, ended, rna_seq, process, location = line.strip().split( "," )
        if ' ' in condition:
            condition = condition.replace( " ", "_" )
        if process == "0":continue
        if rna_seq == "1" and condition not in all_samples:
            all_samples[condition] = {}
        if rna_seq == "0" and condition not in small_rna_samples:
            # print(line)
            small_rna_samples[condition] = {}
        if rna_seq == "1":

            all_samples[condition][Run] = {"bioproject":BioProject,
                                     "condition":condition,
                                     "Date":Date,
                                     "Ended":ended,
                                     "desc":desc,
                                     "read_length":read_length,
                                     "error_corrected":0,
                                     "location_directory":location if location != "-1" else options.raw_data_downloaded_from_NCBI,
                                     "downloaded_from_NCBI":1 if location == "-1" else 0
                                     }
        else:
            small_rna_samples[condition][Run] = {"bioproject":BioProject,
                                     "condition":condition,
                                     "Date":Date,
                                     "Ended":ended,
                                     "desc":desc,
                                     "read_length":read_length,
                                     "error_corrected":0,
                                     "location_directory":location if location != "-1" else options.raw_data_downloaded_from_NCBI,
                                     "downloaded_from_NCBI":1 if location == "-1" else 0
                                     }
    fhr.close()
    options.mrna_md = all_samples
    options.smrna_md = small_rna_samples
    # pprint.pprint(options.mrna_md)
    # pprint.pprint(all_samples)
    with logging_mutex:
        logger_proxy.info( "Metadata information created" )


def readFastaFile( filename ):
    """
    Returns a dictionary with fasta sequences
    """
    if os.path.exists( filename ) == False:return {}
    sequences = {}
    fhr = open( filename, "r" )
    for line in fhr:
        if line[0] == ">":
            sequences[line.strip()[1:].split()[0]] = fhr.readline().strip()
    fhr.close()
    return sequences


def writeFastaFile( filename, sequences ):
    """
    Writes each sequence in the dictionary sequence into the file
    """
    fhw = open( filename, "w" )
    for id in sequences:
        fhw.write( ">" + id + "\n" + sequences[id] + "\n" )
    fhw.close()


def readMultiLineFasta( filename ):
    data = {}
    if os.path.exists( filename ) == False:return {}
    fhr = open( filename, "r" )
    seq, id = "", ""
    line = fhr.readline()
    while True:
        if ">" in line:
            if seq != "":
                data[id] = seq
            id = line.strip()[1:]
            seq = ""
        else:
            seq += line
        line = fhr.readline().strip()
        if not line:break
    data[id] = seq
    fhr.close()
    return data


def readFromRegtoolsOutput( eachinput ):
    SJ_filename_regtools, Run, condition, options = eachinput
    SJ_info = {}
    fhr = open( SJ_filename_regtools, "r" )
    for line in fhr:
        chromosome, chrom_start, chrom_end, name, read_support, strand, thick_start, thick_end, item_rgb, block_count, block_sizes, block_starts = line.strip().split( "\t" )
        if chromosome not in SJ_info:
            SJ_info = {}
        # if chrom_start+"-"+chrom_end not in SJ_info
        SJ_info[chrom_start + "-" + chrom_end] = int( read_support )
    fhr.close()
    return SJ_info, Run, condition


def collectStatsAboutMapping( options ):
    info_about_runs = {}
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:

            if Run not in info_about_runs:
                info_about_runs[Run] = {
                    "total_reads":0,
                                    "Round1":{"uniq_reads_mapped":0, "mm_reads_mapped":0},
                                      "Round2":{"uniq_reads_mapped":0, "mm_reads_mapped":0},
                                      "Round3":{"uniq_reads_mapped":0, "mm_reads_mapped":0},
                                      "round4":{"uniq_reads_mapped":0, "mm_reads_mapped":0}
                                      }
            total_reads, umr, mmr = pullOutMappingInformation( options.output_star + "/" + Run + "_round1_Log.final.out" )
            info_about_runs[Run]["total_reads"] = total_reads
            info_about_runs[Run]["Round1"]["uniq_reads_mapped"] = umr / info_about_runs[Run]["total_reads"]
            info_about_runs[Run]["Round1"]["mm_reads_mapped"] = mmr / info_about_runs[Run]["total_reads"]

            total_reads, umr, mmr = pullOutMappingInformation( options.output_star + "/" + Run + "_round2_Log.final.out" )
            info_about_runs[Run]["Round2"]["uniq_reads_mapped"] = ( info_about_runs[Run]["Round1"]["uniq_reads_mapped"] * info_about_runs[Run]["total_reads"] + umr ) / info_about_runs[Run]["total_reads"]
            info_about_runs[Run]["Round2"]["mm_reads_mapped"] = ( info_about_runs[Run]["Round1"]["mm_reads_mapped"] * info_about_runs[Run]["total_reads"] + mmr ) / info_about_runs[Run]["total_reads"]

            total_reads, umr, mmr = pullOutMappingInformation( options.output_star + "/" + Run + "_round3_Log.final.out" )
            info_about_runs[Run]["Round3"]["uniq_reads_mapped"] = ( info_about_runs[Run]["Round2"]["uniq_reads_mapped"] * info_about_runs[Run]["total_reads"] + umr ) / info_about_runs[Run]["total_reads"]
            info_about_runs[Run]["Round3"]["mm_reads_mapped"] = ( info_about_runs[Run]["Round2"]["mm_reads_mapped"] * info_about_runs[Run]["total_reads"] + mmr ) / info_about_runs[Run]["total_reads"]

            total_reads, umr, mmr = pullOutMappingInformation( options.output_star + "/" + Run + "_round4_Log.final.out" )
            info_about_runs[Run]["round4"]["uniq_reads_mapped"] = ( info_about_runs[Run]["Round3"]["uniq_reads_mapped"] * info_about_runs[Run]["total_reads"] + umr ) / info_about_runs[Run]["total_reads"]
            info_about_runs[Run]["round4"]["mm_reads_mapped"] = ( info_about_runs[Run]["Round3"]["mm_reads_mapped"] * info_about_runs[Run]["total_reads"] + mmr ) / info_about_runs[Run]["total_reads"]

    fhw = open( options.output_star + "/mapping_stats.csv", "w" )
    fhw.write( "Bioproject,Run,condition,Total reads,Round1_umr,Round1_mmr,Round1_tmr,Round2_umr,Round2_mmr,Round2_tmr,Round3_umr,Round3_mmr,Round3_tmr,round4_umr,round4_mmr,round4_tmr\n" )
    for condition in options.mrna_md:
        for Run in options.mrna_md[condition]:
            fhw.write( ",".join( list( map( str, [options.mrna_md[condition][Run]["bioproject"], Run, condition,
             info_about_runs[Run]["total_reads"],
             round( info_about_runs[Run]["Round1"]["uniq_reads_mapped"], 2 ),
             round( info_about_runs[Run]["Round1"]["mm_reads_mapped"], 2 ),
             round( info_about_runs[Run]["Round1"]["uniq_reads_mapped"] + info_about_runs[Run]["Round1"]["mm_reads_mapped"], 2 ),

             round( info_about_runs[Run]["Round2"]["uniq_reads_mapped"], 2 ),
             round( info_about_runs[Run]["Round2"]["mm_reads_mapped"], 2 ),
             round( info_about_runs[Run]["Round2"]["uniq_reads_mapped"] + info_about_runs[Run]["Round2"]["mm_reads_mapped"], 2 ),

             round( info_about_runs[Run]["Round3"]["uniq_reads_mapped"], 2 ),
             round( info_about_runs[Run]["Round3"]["mm_reads_mapped"], 2 ),
             round( info_about_runs[Run]["Round3"]["uniq_reads_mapped"] + info_about_runs[Run]["Round3"]["mm_reads_mapped"], 2 ),

             round( info_about_runs[Run]["round4"]["uniq_reads_mapped"], 2 ),
             round( info_about_runs[Run]["round4"]["mm_reads_mapped"], 2 ),
             round( info_about_runs[Run]["round4"]["uniq_reads_mapped"] + info_about_runs[Run]["round4"]["mm_reads_mapped"], 2 ),
             ] ) ) ) + "\n" )

    fhw.close()


def pullOutMappingInformation( filename ):
    if os.path.exists( filename ) == False:
        return 1, 0, 0
    fhr = open( filename, "r" )
    for line in fhr:
        if "Number of input reads" in line:
            total_reads = line.strip().split( "|" )[-1].strip()
        elif "Uniquely mapped reads number" in line:
            umr = line.strip().split( "|" )[-1].strip()
        elif "Number of reads mapped to multiple loci" in line:
            mmr = line.strip().split( "|" )[-1].strip()
    fhr.close()
    return list( map( int, [total_reads, umr, mmr] ) )


def readAllTranscriptsFromGTFFileInParallel( eachinput ):
    gtf_filename, Run, condition = eachinput
    whole_annotations = {}
    fhr = open( gtf_filename, "r" )
    for line in fhr:
        if line[0] == "#":continue
        chromosome, annotator, structure, start, end, useless1, direction, useless2, desc = line.strip().split( "\t" )
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
        if structure == "transcript": continue
        if transcript_id not in whole_annotations:
            whole_annotations[transcript_id] = {"exons":[],
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
                                                        "chromosome":chromosome,
                                                        "annotator":annotator
                                                        }
        if structure == "exon":
            whole_annotations[transcript_id]["exons"].append( [int( start ), int( end )] )
        elif structure == "CDS":
            whole_annotations[transcript_id]["cds"].append( [int( start ), int( end )] )

        # Construct the introns for each transcript
        if len( whole_annotations[transcript_id]["exons"] ) > 1:
            exon1 = int( whole_annotations[transcript_id]["exons"][-2][1] )
            exon2 = int( whole_annotations[transcript_id]["exons"][-1][0] )
            whole_annotations[transcript_id]["introns"].append( [exon1 + 1, exon2 - 1] )
    fhr.close()

    for transcript_id in whole_annotations:
        whole_annotations[transcript_id]["transcript_start"] = whole_annotations[transcript_id]["exons"][0][0]
        whole_annotations[transcript_id]["transcript_end"] = whole_annotations[transcript_id]["exons"][-1][1]
    for transcript_id in whole_annotations:
        if whole_annotations[transcript_id]["direction"] == "+":
            whole_annotations[transcript_id]["cds_frame"].append( 0 )
            for cds_num, cds in enumerate( whole_annotations[transcript_id]["cds"] ):
                if cds_num == 0:continue
                length_of_CDS_sequence = sum( [cds[1] - cds[0] + 1 for cds in whole_annotations[transcript_id]["cds"][:cds_num]] )
                if length_of_CDS_sequence % 3 == 0:
                    whole_annotations[transcript_id]["cds_frame"].append( 0 )
                elif length_of_CDS_sequence % 3 == 1:
                    whole_annotations[transcript_id]["cds_frame"].append( 2 )
                elif length_of_CDS_sequence % 3 == 2:
                    whole_annotations[transcript_id]["cds_frame"].append( 1 )
        elif whole_annotations[transcript_id]["direction"] == "-":
            whole_annotations[transcript_id]["cds_frame"].append( 0 )
            for cds_num, cds in enumerate( whole_annotations[transcript_id]["cds"][::-1] ):
                if cds_num == 0:continue
                length_of_CDS_sequence = sum( [cds[1] - cds[0] + 1 for cds in whole_annotations[transcript_id]["cds"][::-1][:cds_num]] )
                if length_of_CDS_sequence % 3 == 0:
                    whole_annotations[transcript_id]["cds_frame"].append( 0 )
                elif length_of_CDS_sequence % 3 == 1:
                    whole_annotations[transcript_id]["cds_frame"].append( 2 )
                elif length_of_CDS_sequence % 3 == 2:
                    whole_annotations[transcript_id]["cds_frame"].append( 1 )
            rev = whole_annotations[transcript_id]["cds_frame"][::-1]
            whole_annotations[transcript_id]["cds_frame"] = rev
        else:
            pass
    return whole_annotations, Run, condition


def writeTranscriptsToFile( eachinput ):
    """
    """
    transcript_info, outputfilename, write_introns = eachinput
    fhw = open( outputfilename, "w" )

    for transcript_id in transcript_info:
        chromosome = transcript_info[transcript_id]["chromosome"]
        transcript_start = transcript_info[transcript_id]["exons"][0][0]
        transcript_end = transcript_info[transcript_id]["exons"][-1][1]
        gene_id = ".".join( transcript_id.split( "." )[:2] )
        annotator = transcript_info[transcript_id]["annotator"]
        desc = "gene_id \"" + gene_id + "\"; transcript_id \"" + transcript_id + "\"; "
        desc += "FPKM \"" + str( transcript_info[transcript_id]["FPKM"] ) + "\"; "
        desc += "TPM \"" + str( transcript_info[transcript_id]["TPM"] ) + "\"; "
        desc += "cov \"" + str( transcript_info[transcript_id]["cov"] ) + "\"; "
        row = [chromosome,
         annotator,
         "transcript",
         str( transcript_start ),
         str( transcript_end ),
         "1000",
         transcript_info[transcript_id]["direction"],
         ".",
         desc]
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
             desc]
            fhw.write( "\t".join( row ) + "\n" )
        if write_introns == 1:
            for intron in transcript_info[transcript_id]["introns"]:
                intron_start, intron_end = intron[0], intron[1]
                row = [chromosome,
                 annotator,
                 "intron",
                 str( intron_start ),
                 str( intron_end ),
                 "1000",
                 transcript_info[transcript_id]["direction"],
                 ".",
                 desc]
                fhw.write( "\t".join( row ) + "\n" )
        if "cds" in transcript_info[transcript_id]:
            for cds in transcript_info[transcript_id]["cds"]:
                exon_start, exon_end = cds[0], cds[1]
                row = [chromosome,
                 annotator,
                 "CDS",
                 str( exon_start ),
                 str( exon_end ),
                 "1000",
                 transcript_info[transcript_id]["direction"],
                 ".",
                 desc]
                fhw.write( "\t".join( row ) + "\n" )
    fhw.close()


def translate( seq ):

    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    protein = ""
    if len( seq ) % 3 == 0:
        for i in range( 0, len( seq ), 3 ):
            codon = seq[i:i + 3].upper()
            if codon not in table:
                protein += 'M'
            else:
                protein += table[codon]
    return protein


def splitFasta( input_filename, output_file_prefix, number_of_splits ):
    all_sequences = readFastaFile( input_filename )
    all_sequence_ids = list( all_sequences.keys() )
    all_sequence_ids_split = divide_chunks( all_sequence_ids, len( all_sequence_ids ) // number_of_splits + 1 )
    i = 0
    for i, each_split in enumerate( all_sequence_ids_split ):
        output_filename = f"{output_file_prefix}_{i}.fasta"
        proteins_in_split = {}
        for id in each_split:
            proteins_in_split[id] = all_sequences[id]
        writeFastaFile( output_filename, proteins_in_split )
    return i
