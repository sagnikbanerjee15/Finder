#! /usr/bin/env python
import argparse
import os
import sqlite3
import sys


def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="verify_inputs_to_finder.py",description="Verifies whether all the data are transcriptomic and from the organism under consideration")
    
    # Mandatory arguments
    parser.add_argument("--metadatafile","-mf",help="Please enter the name of the metadata file. Enter 0 in the last column of those samples which you wish to skip processing. The columns should represent the following in order --> BioProject,Run,tissue_group,tissue,description,Date,read_length,ended (PE or SE),directorypath,download,skip. If the sample is skipped it will not be downloaded. Leave the directory path blank if you are downloading the samples. In the end of the run the program will output a csv file with the directory path filled out. Please check the provided csv file for more information on how to configure the metadata file. ",required=True)
    parser.add_argument("--srametadb","-m",help="Enter the location of the SRAmetadb file.",required=True) 
    parser.add_argument("--taxon_id","-t",help="Enter the taxonomic id of the organism. Enter -1 if you are working on a non-model organism or a sub-species for which no taxonomic id exists.",required=True)
    return parser.parse_args()

def readMetaDataFile(options,logger_proxy,logging_mutex):
    all_samples={}
    small_rna_samples={}
    fhr=open(options.metadatafile,"r")
    for line in fhr:
        if "BioProject" in line:continue
        BioProject,Run,condition,desc,Date,read_length,ended,rna_seq,process,location=line.strip().split(",")
        if ' ' in condition:
            condition=condition.replace(" ","_")
        if process=="0":continue
        if rna_seq=="1" and condition not in all_samples:
            all_samples[condition]={}
        if rna_seq=="0" and condition not in small_rna_samples:
            #print(line)
            small_rna_samples[condition]={}
        if rna_seq=="1":
            all_samples[condition][Run]={"bioproject":BioProject,
                                     "condition":condition,
                                     "Date":Date,
                                     "Ended":ended,
                                     "desc":desc,
                                     "read_length":read_length,
                                     "error_corrected":0,
                                     "location_directory":options.raw_data_downloaded_from_NCBI if location=="" else location,
                                     "downloaded_from_NCBI":1 if location=="" else 0
                                     }
        else:
            small_rna_samples[condition][Run]={"bioproject":BioProject,
                                     "condition":condition,
                                     "Date":Date,
                                     "Ended":ended,
                                     "desc":desc,
                                     "read_length":read_length,
                                     "error_corrected":0,
                                     "location_directory":options.raw_data_downloaded_from_NCBI if location=="" else location,
                                     "downloaded_from_NCBI":1 if location=="" else 0
                                     }
    fhr.close()
    options.mrna_md=all_samples
    options.smrna_md=small_rna_samples
    #pprint.pprint(options.mrna_md)
    #pprint.pprint(all_samples)
    with logging_mutex:
        logger_proxy.info("Metadata information created")


def verifyInputs(options):
    """
    """
    if os.path.exists(options.srametadb)==False:
        print("The SRAmetadb does not exists. Program is exiting ")
        sys.exit()
        
    conn = sqlite3.connect(options.srametadb)
    c = conn.cursor()
    study_to_runs={}
    
    perfect=[]
    unknown_source=[]
    wrong_source=[]
    unknown_taxon=[]
    wrong_taxon=[]
    unknown_source_unknown_taxon=[]
    notfound=[]
    for condition in options.mrna_md:
        for run in options.mrna_md[condition]:
            query="""SELECT sra.run_accession,sra.study_accession,sra.library_source,sra.taxon_id FROM sra WHERE run_accession = '"""+run+"""' """
            #print(query)
            flag=0
            for row in c.execute(query):
                run_accession,study_accession,source,taxon_id=row
                taxon_id=str(taxon_id)
                if row==None:
                    flag=1
                    break
                if source is not None and source.upper()=='TRANSCRIPTOMIC' and taxon_id==options.taxon_id:
                    #print(run,"perfect")
                    perfect.append(run)
                elif (source is None or source.upper()!='TRANSCRIPTOMIC') and taxon_id==options.taxon_id:
                    #print(run,"")
                    if source == None:
                        unknown_source.append(run)
                    else:
                        wrong_source.append(run+" "+source.upper())
                elif source is not None and source.upper()=='TRANSCRIPTOMIC' and taxon_id!=options.taxon_id:
                    if taxon_id == "None":
                        unknown_taxon.append(run)
                    else:
                        wrong_taxon.append(run+" "+taxon_id)
                else:
                    unknown_source_unknown_taxon.append(run)
                
            if flag==0:
                notfound.append(run)
    
    print("Perfect runs:\n","\n".join(list(set(perfect))))
    print("Unknown source Correct taxon:\n","\n".join(list(set(unknown_source))))
    print("Wrong source Correct taxon:\n","\n".join(list(set(wrong_source))))
    print("Transcriptomic Unknown taxon:\n","\n".join(list(set(unknown_taxon))))
    print("Transcriptomic Wrong  taxon:\n","\n".join(list(set(wrong_taxon))))
    print("Unknown source Unknown taxon:\n","\n".join(list(set(unknown_source_unknown_taxon))))
    
def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    
    readMetaDataFile(options)
    verifyInputs(options)

if __name__ == "__main__":
    main()
