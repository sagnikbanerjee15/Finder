#ifndef JIGSAWALN_H
#define JIGSAWALN_H

#include <list>
#include <stdint.h>

#include "bwt.h"
#include "splicesitemap.h"

#define BWA_TYPE_NO_MATCH 0
#define BWA_TYPE_UNIQUE 1
#define BWA_TYPE_REPEAT 2
#define BWA_TYPE_MATESW 3

#define SAM_FPD   1 // paired
#define SAM_FPP   2 // properly paired
#define SAM_FSU   4 // self-unmapped
#define SAM_FMU   8 // mate-unmapped
#define SAM_FSR  16 // self on the reverse strand
#define SAM_FMR  32 // mate on the reverse strand
#define SAM_FR1  64 // this is read one
#define SAM_FR2 128 // this is read two
#define SAM_FSC 256 // secondary alignment

#define BWA_AVG_ERR 0.02
#define BWA_MIN_RDLEN 35 // for read trimming

#ifndef bns_pac
#define bns_pac(pac, k) ((pac)[(k)>>2] >> ((~(k)&3)<<1) & 3)
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.1.7"
#endif

using namespace std;

typedef struct {
	bwtint_t w;
	int bid;
} bwt_width_t;

typedef struct {
	uint32_t n_mm:8, n_gapo_t:8, n_gape_t:8, n_gapo_q:8, n_gape_q:8, a:1;
	bwtint_t k, l;
	int score;
} bwt_aln1_t;


/*a hit of a word in the reference*/
typedef struct {
	int wid; /*the id of the word in the query*/
	int64_t pos_q, pos_t; /*start position in the query and target sequence*/
	uint32_t strand:1;  /*strand: 0-forward, 1-reverse complementary */
} jigsaw_word_hit_t;


/*exon-a clump of hits*/

typedef struct {
	uint32_t colinear; //whether hits in the clump are colinear
	uint32_t n_mm:8, n_gapo_t:8, n_gape_t:8, n_gapo_q:8, n_gape_q:8;
	float coverage; //coverage;

	uint32_t strand:1;
	int64_t start_q, start_t, end_q, end_t;
	uint32_t is_first:1, is_last:1; //first or last exon
	//int n_hits;
	list<jigsaw_word_hit_t*> *hits; //shallow copy of hits

	uint32_t exon_id; //for debugging
	//float score;
} jigsaw_exon_t;


typedef struct {
	jigsaw_exon_t *uexon, *dexon; //shallow copy of exons
	//NOTE: shallow copy, don't release memory here

	uint32_t n_mm:8, n_gapo_t:8, n_gape_t:8, n_gapo_q:8, n_gape_q:8;
	//record the number of mismatches near exon junctions
	//since there is no mismatch near hits, these should not be redundant with mismatches
	//recorded in exons
	uint32_t an_mm:8;
	//this is to adjust the mismatches when extending into the exons
	//adjust number of mismatch 

	uint32_t sense_strand:2;
	// strand of the splice site, which might be different from the strand of exons

	int64_t start_t, end_t, start_q, end_q;
	//the last/first nucleotide of the upstream/downstream exon

	double logistic_prob;
	//the logit score for the junction
} jigsaw_junction_t;


typedef struct {

	uint32_t n_mm:8, n_gapo_t:8, n_gape_t:8, n_gapo_q:8, n_gape_q:8, diff;

	uint32_t strand:1; //the strand of the query sequence
	uint32_t sense_strand:2; //the strand of the transcript, as determined by the splice sites
	//strand and sense_strand might be different depending on whether the RNA-seq protocol have strand information

	int64_t start_q, start_t, end_q, end_t;
	double logistic_prob;

	list<jigsaw_junction_t*> *junctions; //shallow copy of junctions
} jigsaw_spliced_aln_t;


typedef struct __jigsaw_spliced_aln_cluster_t {
	//uint32_t strand;
	//int64_t start_q, start_t, end_q, end_t;
	float score;

	list<jigsaw_word_hit_t*>* hits; //shallow copy of hits

	list<jigsaw_exon_t*> *exons;

	list<jigsaw_junction_t*> *junctions;

	uint32_t cluster_id; // for debug

} jigsaw_spliced_aln_cluster_t;



//typedef uint16_t bwa_cigar_t;
/* rgoya: If changing order of bytes, beware of operations like:
 *     s->cigar[0] += s->full_len - s->len;
 */

//#define CIGAR_OP_SHIFT 14
//#define CIGAR_LN_MASK 0x3fff

//since we have introns now, we need more bits for CIGAR
typedef uint32_t jigsaw_cigar_t;

#define CIGAR_OP_SHIFT 28
#define CIGAR_LN_MASK 0xfffffff

#define __cigar_op(__cigar) ((__cigar)>>CIGAR_OP_SHIFT)
#define __cigar_len(__cigar) ((__cigar)&CIGAR_LN_MASK)
#define __cigar_create(__op, __len) (((jigsaw_cigar_t)__op)<<CIGAR_OP_SHIFT | (__len))

#define CIGAR_OP_M FROM_M	//0
#define CIGAR_OP_I FROM_I	//1
#define CIGAR_OP_D FROM_D	//2
#define CIGAR_OP_S FROM_S	//3
#define CIGAR_OP_N	4		//intron


typedef struct {
	uint32_t pos;
	uint32_t n_cigar:15, gap_t:8, gap_q:8, mm:8, strand:1, sense_strand:2, nm:12;
	jigsaw_cigar_t *cigar;
} bwt_multi1_t;



//this is a recursive data structure that represents both a query read and a word in the query read
typedef struct __jigsaw_seq_t {

	/*common variables of the two levels*/

	// for multi-threading only
	int tid;

	char *name;
	ubyte_t *seq, *rseq, *qual;
	uint32_t len:20;
	int clip_len;

	// alignments in SA coordinates
	int n_aln;
	bwt_aln1_t *aln;


	/*variables for bwa_seq_t*/
	// alignment summary
	uint64_t c1:28, c2:28, seQ:8; // number of top1 and top2 hits; single-end mapQ
	// NM and MD tags
	uint32_t full_len:20, nm:12;
	char *md;

	//parameters for the main hit
	bwtint_t sa, pos;
	uint32_t strand:1, type:2, dummy:1, extra_flag:8;
	uint32_t sense_strand:2;
	uint32_t n_mm:8, n_gapo_t:8, n_gape_t:8, n_gapo_q:8, n_gape_q:8, mapQ:8;
	int score;
	int n_cigar;
	jigsaw_cigar_t *cigar;

	// multiple hits
	int n_multi;
	bwt_multi1_t *multi;

	//word information
	int n_words;
	int word_size;
	struct __jigsaw_seq_t *words;

	/*variables for jigsaw_word_t*/
	int wid; //, rwid; //word id
	int offset, roffset;  //offset of the word in the query sequence
	int n_occ; //total occurrence of the word
} bwa_seq_t;


typedef struct __jigsaw_seq_t jigsaw_word_t;
//in case we will define separate data structure later

typedef struct __jigsaw_seq_t jigsaw_anchor_seq_t;
//this is for the single anchor search

#define BWA_MODE_GAPE       0x01
#define BWA_MODE_COMPREAD   0x02
#define BWA_MODE_LOGGAP     0x04
#define BWA_MODE_NONSTOP    0x10
#define BWA_MODE_BAM        0x20
#define BWA_MODE_BAM_SE     0x40
#define BWA_MODE_BAM_READ1  0x80
#define BWA_MODE_BAM_READ2  0x100

#define STRAND_MODE_FORWARD 0x01
#define STRAND_MODE_REVERSE 0X02

typedef struct {
	int word_size;
	int word_max_overlap;
	int max_word_diff;
	int max_word_occ;
	int max_intron_size, min_intron_size;
	int min_anchor, known_junc_min_anchor;
	int min_exon_size;
	int max_overhang;
	int s_mm, s_gapo, s_gape;
	int mode;
	int indel_end_skip, max_del_occ, max_entries, n_batch;
	int max_report_multi;
	float fnr;
	uint32_t max_diff, max_gapo, max_gape;
//	int max_seed_diff, seed_len;
	int n_threads;
	int max_top2;
	int trim_qual;
	char *rg, *junction_file, *regression_file;
	splice_site_map_t *splice_site_map;
	uint32_t verbose;
	uint32_t single_anchor_search;
	uint32_t allow_rep_anchor;
	uint32_t non_denovo_search;
	uint32_t report_best_only;
	double min_logistic_prob;
	int strand_mode;
} gap_opt_t;

#define BWA_PET_STD   1
#define BWA_PET_SOLID 2

typedef struct {
	int max_isize, force_isize;
	int max_occ;
	int n_multi, N_multi;
	int type, is_sw, is_preload;
	double ap_prior;
} pe_opt_t;

struct __bwa_seqio_t;
typedef struct __bwa_seqio_t bwa_seqio_t;

#ifdef __cplusplus
extern "C" {
#endif

	gap_opt_t *gap_init_opt();
	void jigsaw_aln_core(const char *prefix, const char *fn_fa, gap_opt_t *opt, int argc, char *argv[]);
	int bwa_cal_maxdiff(int l, double err, double thres);


	//bwa_seqio_t *bwa_seq_open(const char *fn);
	//bwa_seqio_t *bwa_bam_open(const char *fn, int which);
	//void bwa_seq_close(bwa_seqio_t *bs);
	//void seq_reverse(int len, ubyte_t *seq, int is_comp);
	//bwa_seq_t *bwa_read_seq(bwa_seqio_t *seq, int n_needed, int *n, int is_comp, int trim_qual);
	//void jigsaw_free_read_seq(int n_seqs, bwa_seq_t *seqs);


	//extern gap_stack_t;
	//void jigsaw_cal_sa_reg_gap(bwt_t *const bwt[2], bwa_seq_t *seq, gap_stack_t *stack, const gap_opt_t *opt);

	//void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);

	//in cs2nt.c
	void bwa_cs2nt_core(bwa_seq_t *p, bwtint_t l_pac, const ubyte_t *pac);


	/* rgoya: Temporary clone of aln_path2cigar to accomodate for bwa_cigar_t,
	__cigar_op and __cigar_len while keeping stdaln stand alone */
#include "stdaln.h"

	jigsaw_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar);


#ifdef __cplusplus
}
#endif

#endif
