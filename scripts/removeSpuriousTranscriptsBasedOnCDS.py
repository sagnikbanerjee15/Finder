
from scripts.fileReadWriteOperations import *


def removeSpuriousTranscriptsBasedOnCDS( options, logger_proxy, logging_mutex ):
    """
    """
    inputfilename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS.gtf"
    all_transcripts_info = readAllTranscriptsFromGTFFileInParallel( [inputfilename, "dummy", "dummy"] )[0]
    if os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_high_conf.gtf" ) == True and os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_low_conf.gtf" ) == True:return
    all_transcripts_info_gene_to_transcripts = {}
    # Restructure data structure to group transcripts of same gene together

    for transcript_id in all_transcripts_info:
        gene_id = all_transcripts_info[transcript_id]['gene_id']
        if gene_id not in all_transcripts_info_gene_to_transcripts:
            all_transcripts_info_gene_to_transcripts[gene_id] = {}
        all_transcripts_info_gene_to_transcripts[gene_id][transcript_id] = all_transcripts_info[transcript_id]

    # Low confidence transcripts - Has much lesser CDS length than max(CDS length for all transcripts of the same gene)
    for gene_id in all_transcripts_info_gene_to_transcripts:
        max_CDS_length = 0
        for transcript_id in all_transcripts_info_gene_to_transcripts[gene_id]:
            cds_length = sum( [row[1] - row[0] + 1 for row in all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["cds"]] )
            if max_CDS_length < cds_length:
                max_CDS_length = cds_length
        for transcript_id in all_transcripts_info_gene_to_transcripts[gene_id]:
            cds_length = sum( [row[1] - row[0] + 1 for row in all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["cds"]] )
            if cds_length <= max_CDS_length * 0.5:
                all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["confidence"] = "low"
                # print("Low confidence - Has much lesser CDS length than max(CDS length for all transcripts of the same gene")
            else:
                all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["confidence"] = "high"

    # Transcripts in repeat regions
    cmd = "gffread "
    cmd += " -w " + inputfilename[:-3] + "fasta "
    cmd += " -g " + options.genome
    cmd += " " + inputfilename
    os.system( cmd )

    cmd = "perl " + " -pe '/^>/ ? print \"\\n\" : chomp' "
    cmd += " " + inputfilename[:-3] + "fasta "
    cmd += "| tail -n +2  > " + inputfilename[:-3] + "fasta.temp "
    os.system( cmd )

    cmd = "mv " + inputfilename[:-3] + "fasta.temp " + inputfilename[:-3] + "fasta "
    os.system( cmd )

    transcript_fasta = readFastaFile( inputfilename[:-3] + "fasta" )
    for transcript_id in transcript_fasta:
        if sum( [1 for c in transcript_fasta[transcript_id] if c.upper() == True] ) / len( transcript_fasta[transcript_id] ) < 0.25:
            gene_id = ".".join( transcript_id.split( "." )[:-1] )
            all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["confidence"] = "repeat_region"
            # print("Repeat region")

    # covsplit transcripts without CDS - mark them as low confidence
    for gene_id in all_transcripts_info_gene_to_transcripts:
        for transcript_id in all_transcripts_info_gene_to_transcripts[gene_id]:
            if ( "covsplit" in transcript_id or "isplit" in transcript_id ) and len( all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["cds"] ) == 0:
                all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["confidence"] = "low"
            else:
                all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["confidence"] = "high"
                # print("Low confidence - covsplit")

    # Divide low confident and high confident sets
    high_confidence_genes = {}
    low_confidence_genes = {}
    total_high_confidence_transcripts = 0
    total_low_confidence_transcripts = 0
    for gene_id in all_transcripts_info_gene_to_transcripts:
        for transcript_id in all_transcripts_info_gene_to_transcripts[gene_id]:
            chromosome = all_transcripts_info[transcript_id]["chromosome"]
            if all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["confidence"] == "high":
                high_confidence_genes[transcript_id] = all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]
                total_high_confidence_transcripts += 1
            elif all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["confidence"] == "low" or all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]["confidence"] == "repeat_region":
                low_confidence_genes[transcript_id] = all_transcripts_info_gene_to_transcripts[gene_id][transcript_id]
                total_low_confidence_transcripts += 1

    # print(total_high_confidence_transcripts,total_low_confidence_transcripts)
    writeTranscriptsToFile( [high_confidence_genes, options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_high_conf.gtf", 0] )
    writeTranscriptsToFile( [low_confidence_genes, options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_with_CDS_low_conf.gtf", 0] )

