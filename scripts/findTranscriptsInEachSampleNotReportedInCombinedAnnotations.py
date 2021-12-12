

######################################################################################################################################################################
# File consists of functions for collecting transcripts reported in each of the RNA-Seq Samples
######################################################################################################################################################################

from scripts.alignReads import *
from scripts.fileReadWriteOperations import *
import collections
import multiprocessing


def findTranscriptsInEachSampleNotReportedInCombinedAnnotations( options, logger_proxy, logging_mutex ):
    combined_gtf_filename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined.gtf"
    output_gtf_filename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_missed_transcripts_added.gtf"
    if os.path.exists( output_gtf_filename ) == True:return
    combined_transcript_info = readAllTranscriptsFromGTFFileInParallel( [combined_gtf_filename, "combined", "combined"] )[0]
    number_of_old_transcripts = len( combined_transcript_info )
    combined_transcripts = [transcript_id for transcript_id in combined_transcript_info]
    condition_to_samples_to_unique_transcripts = {}
    pool = multiprocessing.Pool( processes = int( options.cpu ) )
    allinputs = []
    for condition in options.mrna_md:
        condition_to_samples_to_unique_transcripts[condition] = {}
        for Run in options.mrna_md[condition]:
            gtf_filename = options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/psiclass_output_sample_" + Run + ".gtf"
            allinputs.append( [gtf_filename, Run, condition] )
            condition_to_samples_to_unique_transcripts[condition][Run] = {"transcript_info":[], "transcript":[]}
    results = pool.map( readAllTranscriptsFromGTFFileInParallel, allinputs )
    retain_these_transcripts = {}
    for condition in options.mrna_md:
        retain_these_transcripts[condition] = []
        for Run in options.mrna_md[condition]:
            condition_to_samples_to_unique_transcripts[condition][Run]["transcript_info"] = [result[0] for result in results if result[1] == Run and result[2] == condition][0]
            transcript_info = [result[0] for result in results if result[1] == Run and result[2] == condition][0]
            condition_to_samples_to_unique_transcripts[condition][Run]["transcript"] = [transcript_id for transcript_id in transcript_info if len( transcript_info[transcript_id]["exons"] ) > 1]
            condition_to_samples_to_unique_transcripts[condition][Run]["transcript"] = set( condition_to_samples_to_unique_transcripts[condition][Run]["transcript"] ) - set( combined_transcripts )
            # print(condition_to_samples_to_unique_transcripts[condition][Run]["transcript"][:5])

        transcripts_present_in_condition = []
        for Run in options.mrna_md[condition]:
            transcripts_present_in_condition.extend( condition_to_samples_to_unique_transcripts[condition][Run]["transcript"] )
        counter = collections.Counter( transcripts_present_in_condition )
        for transcript_id in counter:
            if counter[transcript_id] == len( options.mrna_md[condition] ):
                retain_these_transcripts[condition].append( transcript_id )

    # Adding the retained transcripts into the combined annotation
    for condition in retain_these_transcripts:
        for transcript_id in retain_these_transcripts[condition]:
            Run = list( options.mrna_md[condition].keys() )[0]
            combined_transcript_info[transcript_id] = condition_to_samples_to_unique_transcripts[condition][Run]["transcript_info"][transcript_id]

    writeTranscriptsToFile( [combined_transcript_info, output_gtf_filename, 0] )
    number_of_new_transcripts = len( combined_transcript_info )
    with logging_mutex:
        logger_proxy.info( f" {number_of_new_transcripts-number_of_old_transcripts} new transcripts added after processing each RNA-Seq sample" )
