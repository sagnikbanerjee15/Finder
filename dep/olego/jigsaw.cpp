#include <stdio.h>
//#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
//#include "bwtaln.h"
#include "jigsawaln.h"
#include "bwtgap.h"
#include "utils.h"


int main (int argc, char *argv[])
{
	int c, opte = -1;
	gap_opt_t *opt;
	//char *bwa_rg = 0;
	opt = gap_init_opt();

	/*
	static struct option long_options [] = {
		{"word-size", 		required_argument, 0, 'w'},
		{"max-word-occ",	required_argument, 0, 'W'},
		{"junction-file",	required_argument, 0, 'j'},
		{"max-diff", 		required_argument, 0, 'n'},
		{"max-gapo", 		required_argument, 0, 'o'},
		{"max-gape", 		required_argument, 0, 'e'},
		{"max-intron-size",	required_argument, 0, 'I'},
		{"indel-end-skip", 	required_argument, 0, 'i'},
		{"gape-max-occ", 	required_argument, 0, 'd'},
		{"max-entries", 	required_argument, 0, 'm'},
		{"num-threads",		required_argument, 0, 't'},
		{"penalty-mismatch",required_argument, 0, 'M'},
		{"penalty-gapo",	required_argument, 0, 'O'},
		{"penalty-gape",	required_argument, 0, 'E'},
		{"repeat",			required_argument, 0, 'R'},
		{"rg-line",			required_argument, 0, 'r'},
		{"output-file",     required_argument, 0, 'f'},
		{"color-space",		no_argument, 0, 'c'},
		{"log-gap",			no_argument, 0, 'L'},
		{"none-stop",		no_argument, 0, 'N'},
		{"bam",				no_argument, 0, 'b'},
		{"single",			no_argument, 0, '0'},
		{"read1",			no_argument, 0, '1'},
		{"read2",			no_argument, 0, '2'},
		{"verbose",			no_argument, 0, 'v'},
		{0, 0, 0, 0}
	};
*/
	static struct option long_options [] = {

			/*input options*/
			{"color-space",		no_argument, 		0, 'c'},
			{"bam",				no_argument, 		0, 0},
			{"single",			no_argument, 		0, '0'},
			{"read1",			no_argument,		0, '1'},
			{"read2",			no_argument, 		0, '2'},
			{"rg-line",			required_argument, 	0, 0},

			/*basic options*/
			{"output-file",     required_argument, 	0, 'o'},
			{"word-size", 		required_argument, 	0, 'w'},
			{"max-word-occ",	required_argument, 	0, 'W'},
			{"max-total-diff", 	required_argument, 	0, 'M'},
			{"max-word-diff",	required_argument, 	0, 'm'},
			{"max-intron",		required_argument, 	0, 'I'},
			{"min-intron",		required_argument, 	0, 'i'},
			{"min-exon", 		required_argument,  0, 'e'},
			{"min-anchor",		required_argument, 	0, 'a'},
			{"known-junc-min-anchor",   required_argument,      0, 'k'},
			{"regression-model",		required_argument,      0, 'r'},
			{"junction-file",	required_argument, 	0, 'j'},
			{"non-denovo",      no_argument,        0, 'n'},
			{"num-threads",		required_argument, 	0, 't'},
//			{"best",		no_argument,		0, 'b'},
			{"verbose",			no_argument, 		0, 'v'},

			/*advanced options*/
			{"non-single-anchor",	no_argument,		0, 0},
			{"allow-rep-anchor",   no_argument,            0, 0},
			{"strand-mode",         required_argument,      0, 0},
			{"word-max-overlap",	required_argument,	0, 0},
			{"max-multi",           required_argument,      0, 0},
			{"min-logistic-prob",         required_argument,      0, 0},
			{"max-overhang",    required_argument,  0, 0},
			{"max-gapo", 		required_argument, 	0, 0},
			{"max-gape", 		required_argument, 	0, 0},
			{"indel-end-skip", 	required_argument, 	0, 0},
			{"gape-max-occ", 	required_argument, 	0, 0},
			{"penalty-mismatch",required_argument, 	0, 0},
			{"penalty-gapo",	required_argument, 	0, 0},
			{"penalty-gape",	required_argument, 	0, 0},
			{"log-gap",			no_argument, 		0, 'L'},
			{"num-reads-batch",         required_argument,      0, 0},
			{"max-entries", 	required_argument, 	0, 0},
			{"repeat",			required_argument, 	0, 0},
			{"none-stop",		no_argument, 		0, 0},

			{0, 0, 0, 0}
		};


	int option_index = 0;

	while ((c = getopt_long (argc, argv,
				"c012o:w:W:M:m:I:i:e:a:j:r:t:k:vLsn",
				long_options, &option_index)) >= 0) {
		switch (c) {
		case -1 : break; //the end of the options

		case 'c': opt->mode &= ~BWA_MODE_COMPREAD; break;
		case '0': opt->mode |= BWA_MODE_BAM_SE; break;
		case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
		case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
		case 'o': freopen(optarg, "w", stdout); break;
		case 'w': opt->word_size = atoi(optarg); break;
		case 'W': opt->max_word_occ = atoi(optarg); break;
		case 'M':
			if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
			else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			break;
		case 'm': opt->max_word_diff = atoi(optarg); break;
		case 'I': opt->max_intron_size = atoi(optarg); break;
		case 'i': opt->min_intron_size = atoi(optarg); break;
		case 'e': opt->min_exon_size = atoi (optarg); break;
		case 'a': opt->min_anchor = atoi(optarg); break;
		case 'k': opt->known_junc_min_anchor = atoi(optarg); break;
//		case 's': opt->single_anchor_search = 0; break;
		case 'n': opt->non_denovo_search = 1; break;
		case 'r': opt->regression_file = optarg; break;
		case 'j': opt->junction_file = optarg; break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'v': opt->verbose = 1; break;


		case 0  : //long options without short representations
			if (long_options[option_index].flag != 0)
				break;

			if (strcmp (long_options[option_index].name, "bam" ) == 0) {
				opt->mode |= BWA_MODE_BAM;
			}
			else if (strcmp (long_options[option_index].name, "rg-line" ) == 0) {
				opt->rg = optarg;
			}
			else if (strcmp (long_options[option_index].name, "max-overhang" ) == 0) {
				opt->max_overhang = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "max-gapo" ) == 0) {
				opt->max_gapo = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "max-gape" ) == 0) {
				opte = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "indel-end-skip" ) == 0) {
				opt->indel_end_skip = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "gape-max-occ" ) == 0) {
				opt->max_del_occ = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "penalty-mismatch" ) == 0) {
				opt->s_mm = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "penalty-gapo" ) == 0) {
				opt->s_gapo = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "penalty-gape" ) == 0) {
				opt->s_gape = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "log-gap" ) == 0) {
				opt->mode |= BWA_MODE_LOGGAP;
			}
			else if (strcmp (long_options[option_index].name, "num-reads-batch" ) == 0) {
			        opt->n_batch = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "max-multi" ) == 0) {
			                                    opt->max_report_multi = atoi(optarg);
			}
									    							    
			else if (strcmp (long_options[option_index].name, "max-entries" ) == 0) {
				opt->max_entries = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "repeat" ) == 0) {
				opt->max_top2 = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "none-stop" ) == 0) {
				opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff;
			}
			else if (strcmp (long_options[option_index].name, "min-logistic-prob") == 0 ) {
				opt->min_logistic_prob = atof (optarg);
			}
			else if (strcmp (long_options[option_index].name, "strand-mode") == 0) {
				opt->strand_mode = atoi (optarg);
			}
			else if (strcmp (long_options[option_index].name, "non-single-anchor") == 0) {
				opt->single_anchor_search = 0;
			}
			else if (strcmp (long_options[option_index].name, "word-max-overlap") ==0 ) {
				opt->word_max_overlap = atoi (optarg);
			}
			else if (strcmp (long_options[option_index].name, "allow-rep-anchor") == 0) {
				opt->allow_rep_anchor = 1;
			}

			//	int tmp = optarg;
			//}
			break;


		default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\nOLego: version %s\n", PACKAGE_VERSION);
		fprintf(stderr, "Compiled at %s, %s\n", __TIME__, __DATE__);
		fprintf(stderr, "Usage:   olego [options] <prefix> <in.fastx>\n\n");

		fprintf(stderr, "[Arguments]\n");
		fprintf(stderr, "<prefix>	The base name of the reference sequence index.\n");
		fprintf(stderr, "<in.fastx>	Fasta or fastq input file. (gz file allowed, - for STDIN)\n");
		
		/*
		fprintf(stderr, "[\ninput options]\n");
		fprintf(stderr, " -c,--color-space             input sequences are in the color space\n");
		fprintf(stderr, " --bam                        the input read file is in the BAM format\n");
		fprintf(stderr, " -0,--single                  use single-end reads only (effective with -b)\n");
		fprintf(stderr, " -1,--read1                   use the 1st read in a pair (effective with -b)\n");
		fprintf(stderr, " -2,--read2                   use the 2nd read in a pair (effective with -b)\n");
		fprintf(stderr, " -r,--rg-line          TXT    RG line\n");
		*/

		fprintf(stderr, "\n[basic options]\n");
		fprintf(stderr, " -o,--output-file      FILE   Output file [stdout]\n");

		fprintf(stderr, " -j,--junction-file    FILE   BED file for known junctions \n");
		fprintf(stderr, " -n,--non-denovo              Disable de novo junction search \n");
		fprintf(stderr, " -t,--num-threads      INT    Number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, " -r,--regression-model FILE   The file with the logistic regression model\n");
		fprintf(stderr, " -M,--max-total-diff   INT    Max #diff (int)[%d] or missing prob under %.2f err rate (float) [%.2f]\n",opt->max_diff, BWA_AVG_ERR, opt->fnr);
		fprintf(stderr, " -w,--word-size        INT    Seed size in junction search [%d]\n", opt->word_size);
		fprintf(stderr, " -W,--max-word-occ     INT    Max #occurrence of a seed to be used in the further steps  [max (%d * 3^(%d-word_size), 300) ]\n", 1000, opt->word_size);
		fprintf(stderr, " -m,--max-word-diff    INT    Max #diff allowed a seed [%d]\n", opt->max_word_diff);
		fprintf(stderr, " -I,--max-intron       INT    Max intron size for de novo search [%d]\n", opt->max_intron_size);
		fprintf(stderr, " -i,--min-intron       INT    Min intron size for de novo search [%d]\n", opt->min_intron_size);
		fprintf(stderr, " -e,--min-exon         INT    Min exon size [%d]\n", opt->min_exon_size);
		fprintf(stderr, " -a,--min-anchor       INT    Min anchor size in de novo single-anchor junction search [%d]\n", opt->min_anchor);
		fprintf(stderr, " -k,--known-junc-min-anchor   INT Min anchor size for a known junction [%d]\n",opt->known_junc_min_anchor);
//		fprintf(stderr, " -b,--best                    only report the best alignments\n");
		fprintf(stderr, " -v,--verbose                 Verbose mode\n");

		fprintf(stderr, "\n[advanced options]\n");
		fprintf(stderr, " --word-max-overlap    INT    Max # of overlaps between seeds in the seeding step. [%d] \n",opt->word_max_overlap ) ;
		fprintf(stderr, " --non-single-anchor          Disable single-anchor de novo junction search \n");
		fprintf(stderr, " --allow-rep-anchor           Allow anchors with simple repetitive sequences \n");
		fprintf(stderr, " --strand-mode         INT    Strand searching mode, 1:forward, 2: reverse, 3 both. [%d]\n", opt->strand_mode);
		fprintf(stderr, " --max-multi           INT    Max # of multi alignments reported [%d]\n", opt->max_report_multi);
		fprintf(stderr, " --min-logistic-prob   FLOAT  Min logistic probability required for an alignment, in the range of [0,1) [%.2f]\n", opt->min_logistic_prob);
		fprintf(stderr, " --max-overhang        INT    Max # of overhanging nucleotide allowed near junctions [%d]\n", opt->max_overhang);
		fprintf(stderr, " --max-gapo            INT    Max number or fraction of gap opens [%d]\n", opt->max_gapo);
		fprintf(stderr, " --max-gape            INT    Max number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr, " --indel-end-skip      INT    Do not put an indel within INT nt towards the ends [%d]\n", opt->indel_end_skip);
		fprintf(stderr, " --gape-max-occ        INT    Max # occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
		fprintf(stderr, " --penalty-mismatch    INT    Mismatch penalty [%d]\n", opt->s_mm);
		fprintf(stderr, " --penalty-gapo        INT    Gap open penalty [%d]\n", opt->s_gapo);
		fprintf(stderr, " --penalty-gape        INT    Gap extension penalty [%d]\n", opt->s_gape);
		fprintf(stderr, " --log-gap                    Log-scaled gap penalty for long deletions\n");
		fprintf(stderr, " --num-reads-batch     INT    Number of reads per batch [%d]\n", opt->n_batch);
//		fprintf(stderr, " --max-entries         INT    Max entries in the queue [%d]\n", opt->max_entries);
//		fprintf(stderr, " --repeat              INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
		fprintf(stderr, " --none-stop                  Non-iterative mode: search for all n-difference hits (slooow)\n");

		//fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
		//fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		//fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);

		fprintf(stderr, "\n");
		fprintf(stderr, "Please find detailed manual in README or online at http://ngs-olego.sourceforge.net/doc/\n\n");
		return 1;
	}
	if (opt->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
			if (l != k) fprintf(stderr, "[olego_aln] %dnt reads: max_diff = %d\n", i, l);
			k = l;
		}
	}
	fprintf(stderr,"Command used: ");
	for(int i = 0; i < argc; i++ ){
		fprintf(stderr,"%s ", argv[i]);
	}
	fprintf(stderr,"\n");
	if (opt->non_denovo_search && (opt->junction_file == 0)) {
	    fprintf(stderr,"[olego_aln] Warning! Non-denovo search mode without junction annotation, no splice mapping will be reported!\n" );
	}
	//set max_word_occ based on the word size
	if (opt->max_word_occ == -1) { 
		//if user didn't specify this
		opt->max_word_occ =  1000 * pow (3, 14 - opt->word_size );
		if (opt->max_word_occ < 300 ) opt->max_word_occ = 300;
	}
	
	jigsaw_aln_core(argv[optind], argv[optind+1], opt, argc, argv);
	free(opt);
	return 0;
}

