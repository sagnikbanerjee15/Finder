#!/usr/bin/env python3

from inspect import currentframe, getframeinfo
import argparse
import re
import os
import errno
import shutil
import subprocess
from random import choice
from string import ascii_uppercase
import sys

__author__ = "Katharina J. Hoff and Jens Keilwagen"
__copyright__ = "Copyright 2020. All rights reserved."
__license__ = "Artistic Licsense"
__version__ = "1.0.0"
__credits__ = "Maria Hartmann, Ingo Bulla, Mario Stanke"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "production"

parser = argparse.ArgumentParser(
    description='Assembly-free construction of UTRs from short read ' +
                'RNA-Seq data on the basis of coding sequence annotation, ' +
                'produces a gtf file with genes with UTRs and a list of ' +
                'transcript IDs for transcripts that have both 5\'- and ' +
                '3\'-UTR.')
parser.add_argument('-b', '--bam', required=True, type=str, nargs="+",
                    help='File with short RNA-Seq read to genome spliced ' +
                    'alignments in BAM format.')
parser.add_argument('-t', '--gtf', required=True, type=str,
                    help='File with AUGUSTUS CDS-based gene models in GTF ' +
                    'format')
parser.add_argument('-g', '--genome', required=True, type=str,
                    help="Corresponding genome file in FASTA format")
parser.add_argument('-o', '--outfile_name_stem', required=True, type=str,
                    help='Output file name stem for the two files: ' +
                    'AUGUSTUS CDS gene models plus ' +
                    'GeMoMa ' + 'generated UTRs in gtf format and a list ' +
                    'of transcript IDs for those transcripts that have both ' +
                    '5\'- and 3\'-UTR')
parser.add_argument('-d', '--outdir', required=False, type=str, default=".",
                    help='Directory in which a subdirectory for temporary ' +
                    'files of this script run will be created (and ' +
                    'removed). Default ' +
                    'is the current working directory. Must be writable.')
parser.add_argument('-q', '--verbosity', required=False, type=int, default=0,
                    help='Logging verbosity given as positive integer ' +
                    '(default is 0).')
parser.add_argument('-c', '--cores', required=False, type=int, default=1,
                    help='Number of cores for samtools sort processes')
parser.add_argument('-s', '--SAMTOOLS_PATH', required=False, type=str,
                    help="Path to samtools executable")
parser.add_argument('-a', '--AUGUSTUS_SCRIPTS_PATH', required=False, type=str,
                    help="Path to AUGUSTUS scripts")
parser.add_argument('-j', '--JAVA_PATH', required=False, type=str,
                    help="Path to java executable")
parser.add_argument('-m', '--GeMoMaJar', required=False, type=str,
                    help="GeMoMa jar file")
parser.add_argument('-v', '--version', action='version',
                    version='%(prog)s ' + __version__)

args = parser.parse_args()

''' Check whether args.outdir exists and is writable '''

if not os.path.isdir(args.outdir):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': specified argument for --outdir ' +
          "(" + args.outdir + ") is not a directory!")
    exit(1)
elif not os.access(args.outdir, os.W_OK):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': specified argument for --outdir ' +
          "(" + args.outdir + ") is not writable!")
    exit(1)

''' Check whethere GeMoMa_temp exists, this indicates another interfering 
job already running '''

if os.path.isdir('GeMoMa_temp'):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': directory GeMoMa_temp exists, ' +
          'this indicates that an interfering job may currently be ' +
          'running. Either wait until that job has completed, or ' +
          'delete GeMoMa_temp if you are sure that there is no ' +
          'such job.')
    exit(1)


tmp_dir = args.outdir + "/gushr-" + \
    ''.join(choice(ascii_uppercase) for i in range(12)) + "/"
try:
    os.makedirs(tmp_dir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

''' Find bash grep '''

grep_tool = shutil.which('grep')
if grep_tool is None:
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': ' + "Unable to locate bash tool 'grep'")
    print('grep is part of most Linux distributions. On Ubuntu, it is ' +
          'part of the package coreutils. Try re-installing your bash if ' +
          'grep is missing on your system.')
    quit(1)

''' Find java '''

java = None

if args.JAVA_PATH is not None:
    if os.access(args.JAVA_PATH + "/java", os.X_OK):
        java = args.JAVA_PATH + "/java"
else:
    java = shutil.which('java')
    if java is not None:
        if not(os.access(java, os.X_OK)):
            java = None


if java is None:
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': ' +
          "Unable to locate java on your system. ")
    print('Please install java version 1.8.')
    quit(1)

''' Check java version for GeMoMa '''

jv = os.popen(
    "java -version 2>&1 | " + grep_tool + " 'version' 2>&1 | " +
    "awk -F\\\" '{ split($2,a,\".\"); print a[1]\".\"a[2]}'").read()
if not jv == "1.8\n":
    frameinfo = getframeinfo(currentframe())
    print('Warning in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': java version on your system is ' +
          jv + ". This script has only been tested with java version 1.8.")

''' Check whether custom GeMoMa jar is present '''

jar = None
if args.GeMoMaJar:
    jar=args.GeMoMaJar
else:
    this_script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    jar = this_script_path + "/GeMoMa-1.6.2.jar"

if not os.path.isfile(jar):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': custom GeMoMa jar ' +
          jar + ' is missing!')
    exit(1)

''' Find samtools (if bam file provided) '''

samtools = ""

if args.verbosity > 0:
    print("Searching for samtools:")
if args.SAMTOOLS_PATH:
    samtools = args.SAMTOOLS_PATH + "/samtools"
    if not(os.access(samtools, os.X_OK)):
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + samtools + " is not executable!")
        exit(1)
    else:
        if args.verbosity > 0:
            print("Will use " + samtools)
else:
    if shutil.which("samtools") is not None:
        samtools = shutil.which("samtools")
        if args.verbosity > 0:
            print("Will use " + samtools)
    else:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': '
              + "Unable to locate samtools binary!")
        print("samtools is available as package in many Linux distributions.")
        print("For example, on Ubuntu, try installing with:")
        print("\"sudo apt install samtools\"")
        print("If samtools is unavailable as a package, you can obtain it " +
              "from github at:")
        print("https://github.com/samtools/samtools")
        exit(1)

''' Find AUGUSTUS script gtf2gff.pl '''

gtf2gff = ""

if args.verbosity > 0:
    print("Searching for gtf2gff.pl:")
if args.AUGUSTUS_SCRIPTS_PATH:
    gtf2gff = args.AUGUSTUS_SCRIPTS_PATH + "/gtf2gff.pl"
    if not(os.access(gtf2gff, os.X_OK)):
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + gtf2gff + " is not executable!")
        exit(1)
    else:
        if args.verbosity > 0:
            print("Will use " + gtf2gff)
else:
    if shutil.which("gtf2gff.pl") is not None:
        gtf2gff = shutil.which("gtf2gff.pl")
        if args.verbosity > 0:
            print("Will use " + gtf2gff)
    else:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': '
              + "Unable to locate gtf2gff.pl!")
        print("gtf2gff.pl is part of AUGUSTUS scripts.")
        print("You can obtain it " +
              "from github with:")
        print("git clone https://github.com/Gaius-Augustus/Augustus.git")
        print("Compilation and full installation of AUGUSTUS is not " +
              "required for excuting this script. You only need to add " +
              "the missing script to your $PATH.")
        exit(1)

''' Find bash sort '''

sort_tool = shutil.which('sort')
if sort_tool is None:
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': ' + "Unable to locate bash tool 'sort'")
    print('sort is part of most Linux distributions. On Ubuntu, it is ' +
          'part of the package coreutils. Try re-installing your bash if ' +
          'sort is missing on your system.')
    quit(1)

''' ******************* BEGIN FUNCTIONS ************************************'''


''' Function that writes subprocess byte object to flat file '''


def write_byteobj(byte_obj, outfile):
    try:
        with open(outfile, 'w') as byteobj_handle:
            byteobj_handle.write(byte_obj.decode('utf-8'))
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + outfile +
              " for writing!")
        quit(1)


''' Function that runs a simple subprocess, several attempts in case of failure
(because samtools sometimes fails at random'''


def run_simple_process(args_lst):
    max_attempts = 4
    attempts = 0
    while attempts < max_attempts:
        try:
            if args.verbosity > 0:
                print("Trying to execute the following command:")
                print(" ".join(args_lst))
            result = subprocess.run(
                args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if args.verbosity > 0:
                print("Suceeded in executing command.")
            if(result.returncode == 0):
                return(result)
            else:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' +
                      "Return code of subprocess was " +
                      str(result.returncode) + str(result.args))
                attempts += 1
                if attempts == max_attempts:
                    quit(1)
                else:
                    print("Will try again...")
        except subprocess.CalledProcessError as prcsexc:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed executing: ",
                  " ".join(prcsexc.args))
            print("Error code: ", prcsexc.returncode, prcsexc.output)
            attempts += 1
            if attempts == max_attempts:
                quit(1)
            else:
                print("Will try again...")


''' Function that runs a subprocess with input from STDIN '''


def run_process_stdinput(args_lst, byte_obj):
    try:
        if args.verbosity > 0:
            print("Trying to execute the following command with input from " +
                  "STDIN:")
            print(" ".join(args_lst))
        result = subprocess.run(args_lst, input=byte_obj,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if args.verbosity > 0:
            print("Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': '
                  + "run_process_stdinput: return code of subprocess was "
                  + str(result.returncode))
            quit(1)
    except subprocess.CalledProcessError as prcsexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' +
              "Failed executing: ", " ".join(prcsexc.args))
        print("Error code: ", prcsexc.returncode, prcsexc.output)
        quit(1)


''' Function that converts a bam file to a bedgraph '''


def bam_to_bedgraph(bams):
    # currently only unstranded libraries
    # sort bams
    i = 0
    for bam_file in bams:
        bam_sorted_file = tmp_dir + 'rnaseq_' + str(i) + '_s.bam'
        subprcs_args = [samtools, "sort", '-@',
                        str(args.cores), bam_file, "-o", bam_sorted_file]
        run_simple_process(subprcs_args)
        i += 1
    # merge bams
    merged_file = tmp_dir + 'rnaseq_merged.bam'
    subprcs_args = [samtools, "merge", merged_file]
    for j in range(0, i):
        subprcs_args.append(tmp_dir + 'rnaseq_' + str(j) + '_s.bam')
    run_simple_process(subprcs_args)
    # generate bedgraph file
    subprcs_args = [java, '-jar', jar, 'CLI', 'ERE', 'm=' +
                    merged_file, 'u=true', 'c=true', 'outdir='+tmp_dir]
    run_simple_process(subprcs_args)
    return tmp_dir + 'coverage.bedgraph', tmp_dir + 'introns.gff'


''' Function that finds complete genes in a gtf file
(they have both start_codon and stop_codon) '''


def gtf_filter_complete(gtf_file):
    has_start = {}
    has_stop = {}
    try:
        with open(gtf_file, "r") as gtf_handle:
            for line in gtf_handle:
                if re.search(r'transcript_id "([^"]+)";', line):
                    txid = re.search(
                        r'transcript_id "([^"]+)";', line).group(1)
                    if re.search(r'\tstart_codon\t', line):
                        has_start[txid] = 1
                    elif re.search(r'\tstop_codon\t', line):
                        has_stop[txid] = 1
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              gtf_file + " for reading!")
        quit(1)
    has_both = {}  # txids
    for key in has_start.keys():
        if key in has_stop:
            has_both[key] = 1
    has_both_gids = {}
    for key in has_both.keys():
        gid = re.sub(r'([^.]+)\.[^.]+', r'\1', key)
        has_both_gids[gid] = 1
    filtered_file = tmp_dir + 'complete.gtf'
    try:
        with open(filtered_file, "w") as compl_handle:
            try:
                with open(gtf_file, "r") as gtf_handle:
                    for line in gtf_handle:
                        if re.search(r'\tgene\t', line):
                            gid = re.search(r'\t([^\t]+)\n', line).group(1)
                            if gid in has_both_gids:
                                compl_handle.write(line)
                        elif re.search(r'\ttranscript\t', line):
                            if re.search('\ttranscript_id \"[^"]+\"', line):
                                txid = re.search('\ttranscript_id \"([^"]+)\"', line).group(1)
                            else:
                                txid = re.search(r'\t([^\t]+)\n', line).group(1)
                            if txid in has_both:
                                compl_handle.write(line)
                        elif re.search(r'transcript_id "[^"]+";', line):
                            txid = re.search(
                                r'transcript_id "([^"]+)";', line).group(1)
                            if txid in has_both:
                                compl_handle.write(line)
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " +
                      gtf_file + " for reading!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              compl_file + " for writing!")
        quit(1)
    return filtered_file


''' Function that converts AUGUSTUS gtf format to gff3 format '''


def gtf2gff3(gtf_file):
    try:
        with open(gtf_file, "r") as gtf_handle:
            gtf_data = gtf_handle.read()
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              gtf_file + " for reading!")
        quit(1)
    gff3_file = tmp_dir + 'complete.gff3'
    subprcs_args = [gtf2gff, '--out=' + gff3_file, '--gff3']
    run_process_stdinput(subprcs_args, gtf_data.encode('utf-8'))
    print("Done")
    return(gff3_file)


''' Function that makes the AUGUSTUS gff3 compatible with what GeMoMa
expects from gff3 format '''


def augustus_gff3_to_gemoma_gff3(aug_gff3_file):
    gemoma_like_gff3_file = tmp_dir + 'complete_gemoma_like.gff3'
    try:
        with open(aug_gff3_file, 'r') as aug_handle:
            try:
                with open(gemoma_like_gff3_file, 'w') as gemoma_handle:
                    for line in aug_handle:
                        line = line.strip('\n')
                        line_elements = re.split(r'\t', line)
                        if(re.search(r'\tmRNA\t.*ID=([^;]*);', line)):
                            gid = re.search(
                                r'\tmRNA\t.*ID=([^;]*);', line).group(1)
                            nCds = 1
                        if(re.search(r'\tmRNA\t', line)):
                            for i in range(0, 8):
                                if i != 2:
                                    gemoma_handle.write(
                                        line_elements[i] + '\t')
                                else:
                                    gemoma_handle.write("gene\t")
                            gemoma_handle.write(
                                "ID=" + gid + ";transcripts=1;complete=1;" +
                                "maxEvidence=1;maxTie=0.3333\n")
                            for i in range(0, 8):
                                if i != 2:
                                    gemoma_handle.write(
                                        line_elements[i] + '\t')
                                else:
                                    gemoma_handle.write("prediction\t")
                            gemoma_handle.write(
                                "ID=" + gid + "_R1;ref-gene=NA;AA=NA;" +
                                "score=NA;tae=NA;tde=NA;tie=NA;" +
                                "minSplitReads=0;start=M;stop=*;" +
                                "evidence=1;Parent=" + gid + "\n")
                        elif(re.search(r'\tCDS\t', line)):
                            for i in range(0, 8):
                                gemoma_handle.write(line_elements[i] + '\t')
                            gemoma_handle.write(
                                "ID=" + gid + ".CDS" + str(nCds) +
                                ";Parent=" + gid + "_R1;\n")
                            nCds = nCds + 1
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " +
                      aug_gff3_file + " for reading!")
                quit(1)

    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              gemoma_like_gff3_file + " for writing!")
        quit(1)
    return(gemoma_like_gff3_file)


''' Function that adds UTRs to the gff2 file (generated by GeMoMa from
AUGUSTUS CDS predictions and RNA-seq coverage bedgraph) '''


def add_utrs_to_gff3(gff3_file, bed_graph, intron_file):
    gff3_utr_file = tmp_dir + "final_annotation.gff"
    subprcs_args = [java, '-jar', jar, 'CLI', 'AnnotationFinalizer',
                    'u=YES', 'g=' + args.genome, 'a=' + gff3_file,
                    'i=' + intron_file, 'c=UNSTRANDED',
                    'coverage_unstranded=' + bed_graph, 'rename=NO',
                    'outdir=' + tmp_dir]
    run_simple_process(subprcs_args)
    return gff3_utr_file


''' Function that extracts UTR features from gff3 output of
GeMoMa, converts to gtf format '''


def gemoma_gff3_to_gtf(gff3_file):
    gtf_file = tmp_dir + 'utrs.gtf'
    try:
        with open(gff3_file, 'r') as gff3_handle:
            try:
                with open(gtf_file, 'w') as gtf_handle:
                    for line in gff3_handle:
                        if re.search(r'_prime_UTR\t', line):
                            line = line.strip('\n')
                            line_elements = re.split(r'\t', line)
                            gitxids = re.search(
                                r'Parent=([^\.]+)\.([^\.]+)_', line).groups()
                            gid = gitxids[0]
                            txid = gitxids[0] + '.' + gitxids[1]
                            gtf_handle.write(
                                line_elements[0] + '\t' + line_elements[1] +
                                '\t')
                            if re.search(r'three_', line):
                                gtf_handle.write("3'-UTR\t")
                            else:
                                gtf_handle.write("5'-UTR\t")
                            for i in range(3, 8):
                                gtf_handle.write(line_elements[i] + ' \t')
                            gtf_handle.write(
                                'transcript_id "' +
                                txid + '"; gene_id "' + gid + '";\n')
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " +
                      gtf_file + " for writing!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              gff3_file + " for reading!")
        quit(1)
    return gtf_file


''' Function that identifies genes with both 5'- and 3'-UTR,
returns the name of a file with a list with tx ids '''


def find_both_utrs(gtf_file):
    both_utrs_lst = args.outfile_name_stem + '_bothutr.lst'
    three_utr = {}
    five_utr = {}
    try:
        with open(gtf_file, 'r') as gtf_handle:
            for line in gtf_handle:
                if re.search(r'\t5\'-UTR\t', line):
                    tx_id = re.search(
                        r'transcript_id "([^"]+)";', line).group(1)
                    five_utr[tx_id] = 1
                elif re.search(r'\t3\'-UTR\t', line):
                    tx_id = re.search(
                        r'transcript_id "([^"]+)";', line).group(1)
                    three_utr[tx_id] = 1
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              gtf_file + " for reading!")
        quit(1)
    both_utr = {}
    for key in three_utr:
        if key in five_utr:
            both_utr[key] = 1
    try:
        with open(both_utrs_lst, 'w') as both_handle:
            for key in both_utr:
                both_handle.write(key + '\n')
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              both_utrs_lst + " for writing!")
        quit(1)
    return


''' Function that merges the intial AUGUSTUS gtf file with the gene models 
with UTRs from GeMoMa, merging is performed because for generating a genbank 
file for training AUGUSTUS, we want to exclude genes in the neighborhood, 
i.e. avoid having CDS in the flanking region, gene and tx feature coordinates
are re-computed from possible novel UTR features '''


def merge_original_with_utrs(original_gtf, utr_gtf):
    all_gtf = tmp_dir + "all_intermediate.gtf"
    final_gtf = args.outfile_name_stem + ".gtf"
    # merge files
    tmp_gtf = ""
    try:
        with open(original_gtf, 'r') as ori_handle:
            for line in ori_handle:
                if not(re.search(r'\tgene\t', line)) and \
                        not(re.search(r'\ttranscript\t', line)) and \
                        not(re.search(r'\tmRNA\t', line)):
                    tmp_gtf += line
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              original_gtf + " for reading!")
        quit(1)
    try:
        with open(utr_gtf, 'r') as utr_handle:
            for line in utr_handle:
                if not(re.search(r'\tgene\t', line)) and \
                        not(re.search(r'\ttranscript\t', line)) and \
                        not(re.search(r'\tmRNA\t', line)):
                    tmp_gtf += line
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              utr_gtf + " for reading!")
        quit(1)

    subprcs_args = [sort_tool, "-k1,1", "-k4,4n", "-s"]
    result = run_process_stdinput(subprcs_args, tmp_gtf.encode('utf-8'))
    write_byteobj(result.stdout, all_gtf)
    # add missing gene and transcript lines
    genes = {}
    try:
        with open(all_gtf, 'r') as all_handle:
            for line in all_handle:
                if not(re.search(r'^#', line)):
                    txid = re.search(r'transcript_id \"([^"]+)\"', line).group(1)
                    gid = re.search(r'gene_id \"([^"]+)\"', line).group(1)
                    if not gid in genes:
                        genes[gid] = {}
                    if not 'txs' in genes[gid]:
                        genes[gid]['txs'] = {}
                    if not txid in genes[gid]['txs']:
                        genes[gid]['txs'][txid] = {}
                    if not 'lines' in genes[gid]['txs'][txid]:
                        genes[gid]['txs'][txid]['lines'] = []
                    if not re.search(r'\texon\t', line) and not re.search(r'\tintron\t', line):
                        genes[gid]['txs'][txid]['lines'].append(line)
                    fields = re.split(r'\t', line)
                    # find gene boundaries
                    if not 'gene' in genes[gid]:
                        genes[gid]['gene'] = {}
                    if not 'start' in genes[gid]['gene']:
                        genes[gid]['gene']['start'] = fields[3]
                        genes[gid]['gene']['seq'] = fields[0]
                        genes[gid]['gene']['strand'] = fields[6]
                    elif fields[3] < genes[gid]['gene']['start']:
                        genes[gid]['gene']['start'] = fields[3]
                    if not 'end' in genes[gid]['gene']:
                        genes[gid]['gene']['end'] = fields[4]
                    elif fields[4] > genes[gid]['gene']['end']:
                        genes[gid]['gene']['end'] = fields[4]
                    # find tx boundaries
                    if not 'tx' in genes[gid]['txs'][txid]:
                        genes[gid]['txs'][txid]['tx'] = {}
                    if not 'start' in genes[gid]['txs'][txid]['tx']:
                        genes[gid]['txs'][txid]['tx']['start'] = fields[3]
                        genes[gid]['txs'][txid]['tx']['seq'] = fields[0]
                        genes[gid]['txs'][txid]['tx']['strand'] = fields[6]
                    elif fields[3] < genes[gid]['txs'][txid]['tx']['start']:
                        genes[gid]['txs'][txid]['tx']['start'] = fields[3]
                    if not 'end' in genes[gid]['txs'][txid]['tx']:
                        genes[gid]['txs'][txid]['tx']['end'] = fields[4]
                    elif fields[4] > genes[gid]['txs'][txid]['tx']['end']:
                        genes[gid]['txs'][txid]['tx']['end'] = fields[4]
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              all_gtf + " for reading!")
        quit(1)
    try:
        with open(final_gtf, "w") as out_handle:
            for gn in genes.keys():
                out_handle.write(genes[gn]['gene']['seq'] + "\tGUSHR\tgene\t" +
                                 genes[gn]['gene']['start'] +
                                 "\t" + genes[gn]['gene']['end'] + "\t.\t" +
                                 genes[gn]['gene']['strand'] + "\t.\t" + gn +
                                 "\n")
                for txid in genes[gn]['txs'].keys():
                    out_handle.write(genes[gn]['txs'][txid]['tx']['seq'] +
                                     "\tGUSHR\ttranscript\t" +
                                     genes[gn]['txs'][txid]['tx']['start'] +
                                     "\t" +
                                     genes[gn]['txs'][txid]['tx']['end'] +
                                     "\t.\t" +
                                     genes[gn]['txs'][txid]['tx']['strand'] +
                                     "\t.\t" + txid + "\n")
                    for l in genes[gn]['txs'][txid]['lines']:
                        out_handle.write(l)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              final_gtf + " for writing!")
        quit(1)

    print("Done")
    return


''' ******************* END FUNCTIONS *************************************'''


bedgraph, introns = bam_to_bedgraph(args.bam)

complete_genes = gtf_filter_complete(args.gtf)

gff3_genes = gtf2gff3(complete_genes)

gff3_gemoma_like_genes = augustus_gff3_to_gemoma_gff3(gff3_genes)

gff3_with_utrs = add_utrs_to_gff3(gff3_gemoma_like_genes, bedgraph, introns)

utrs_gtf = gemoma_gff3_to_gtf(gff3_with_utrs)

find_both_utrs(utrs_gtf)

merge_original_with_utrs(args.gtf, utrs_gtf)

shutil.rmtree(tmp_dir)
shutil.rmtree('GeMoMa_temp')
