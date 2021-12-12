#! /usr/bin/env python3

import os


def checkForNewData( options ):
    """
    Checks whether any new RNA-Seq sample has been requested for processing
    Deletes braker file to enforce rerun
    Deletes combined_split_transcripts_with_bad_SJ_redundancy_removed.gtf to enforce rerun
    """
    if os.path.exists( options.output_directory + "/samples_processed_in_previous_run" ) == False:
        samples_submitted_in_current_run = [Run for condition in options.mrna_md for Run in options.mrna_md[condition]]
        fhw = open( options.output_directory + "/samples_processed_in_previous_run", "w" )
        fhw.write( "\n".join( samples_submitted_in_current_run ) )
        fhw.close()
        return 0
    else:
        samples_processed_in_previous_run = open( options.output_directory + "/samples_processed_in_previous_run", "r" ).read().split( "\n" )[::-1]
        samples_submitted_in_current_run = [Run for condition in options.mrna_md for Run in options.mrna_md[condition]]
        if len( set( samples_processed_in_previous_run ) & set( samples_submitted_in_current_run ) ) != 0:
            return 1  # There are new RNA-Seq samples deposited by user

        if os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined.gtf" ) == True:
            if os.path.exists( options.output_assemblies_psiclass_terminal_exon_length_modified + "/combined/combined_split_transcripts_with_bad_SJ_redundancy_removed.gtf" ) == True:
                if os.path.exists( options.output_braker + "/braker.gtf" ) == True:
                    return 2  # No new RNA-Seq samples and PsiCLASS, FINDER and BRAKER2 runs have completed successfully
                else:
                    return 3  # No new RNA-Seq samples, PsiCLASS, FINDER done but BRAKER NOT done
            else:
                return 4  # No new RNA-Seq samples, PsiCLASS done but FINDER and BRAKER NOT done
        else:
            return 5  # No new RNA-Seq samples, nothing is done


def determineOptimalStartingPoint( options, logger_proxy, logging_mutex ):
    """
    If options.checkpoint is set to -1 then operations will start from the beginning but will skip steps to regenerate files
    For all other values of options.checkpoint, pre-generated data will be overwritten
    """
    flag = checkForNewData( options )
    # options.checkpoint = 0 # Nothing will be deleted and nothing will be regenerated. Computation will start from the beginning and skip the steps already computed.

    if flag == 2:
        pass  # No new RNA-Seq samples. PsiCLASS, FINDER and BRAKER2 runs have completed successfully. No need to change `options.checkpoint` can go ahead with user request
    elif flag == 3:
        if options.checkpoint > 4:
            options.checkpoint = 4  # No new RNA-Seq samples, PsiCLASS, FINDER done but BRAKER NOT done.
    elif flag == 4:
        if options.checkpoint > 3:
            options.checkpoint = 3  # No new RNA-Seq samples, PsiCLASS done but FINDER and BRAKER NOT done
    elif flag == 5:
        if options.checkpoint > 2:
            options.checkpoint = 2

    if options.checkpoint == 0:
        pass  # Do Nothing
    elif options.checkpoint == 1:  # Align reads to reference genome (Will trigger removal of all alignments and start from beginning)
        os.system( f"rm -rf {options.output_star}/*" )  # Delete contents of alignment folder
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/*" )  # Delete contents of assembly folder
        os.system( f"rm -rf {options.output_braker}/*" )  # Delete contents of BRAKER2 folder
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/braker.gtf" )
    elif options.checkpoint == 2:  # Assemble with PsiCLASS (Will remove all assemblies)
        os.system( f"rm -rf {options.output_star}/*_counts_all_info.pkl" )  # Removing genomic read counts file
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/*" )  # Delete contents of assembly folder
        os.system( f"rm -rf {options.output_braker}/*" )  # Delete contents of BRAKER2 folder
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/braker.gtf" )
    elif options.checkpoint == 3:  # Find genes with FINDER (entails changepoint detection)
        os.system( f"rm -rf {options.output_star}/*_counts_all_info.pkl" )  # Removing genomic read counts file
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/combined_*gtf" )  # Delete contents of assembly folder
        os.system( f"rm -rf {options.output_braker}/*" )  # Delete contents of BRAKER2 folder
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/braker.gtf" )
    elif options.checkpoint == 4:  # Predict genes using BRAKER2
        os.system( f"rm -rf {options.output_braker}/*" )  # Delete contents of BRAKER2 folder
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/combined_with_CDS*.gtf" )
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/braker.gtf" )
    elif options.checkpoint == 5:  # Annotate coding regions
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/combined_with_CDS*.gtf" )
    elif options.checkpoint == 6:  # Merge FINDER annotations with BRAKER2 predictions and protein sequences
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/combined_with_CDS_*.gtf" )
        os.system( f"rm -rf {options.output_assemblies_psiclass_terminal_exon_length_modified}/combined/FINDER_BRAKER_PROT.gtf" )
    with logging_mutex:
        logger_proxy.info( f"Starting FINDER from {options.checkpoint} checkpoint" )

