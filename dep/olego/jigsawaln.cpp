#include <list>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif



#include "bntseq.h"
#include "bwaseqio.h"
//#include "bwtaln.h"
#include "bwase.h"
#include "bwtgap.h"
#include "utils.h"
#include "splicesitemap.h"
#include "jigsawaln.h"
#include "kstring.h"
#include "splicescore.h"

#ifdef HAVE_PTHREAD
#define THREAD_BLOCK_SIZE 1024
#include <pthread.h>
static pthread_mutex_t g_seq_lock = PTHREAD_MUTEX_INITIALIZER;
#endif


gap_opt_t *gap_init_opt()
{
	gap_opt_t *o;
	o = (gap_opt_t*)calloc(1, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	   rate. Voilating this requirement will break pairing! */
	o->word_size = 15;
	o->word_max_overlap = 1;
	o->max_word_diff = 0;
	o->max_word_occ =  -1 ; // will reset this in jigsaw.cpp
	o->min_anchor = 8;
	o->known_junc_min_anchor = 5;
	o->max_overhang = 6;
	o->junction_file = 0;
	o->regression_file = 0;
	o->min_exon_size = 9;
	o->max_intron_size = 500000;
	o->min_intron_size = 20;
	o->single_anchor_search = 1;
	o->allow_rep_anchor = 0;
	o->non_denovo_search = 0;
	o->strand_mode = 3;
	o->min_logistic_prob = 0.5;
	o->splice_site_map = 0;
	o->s_mm = 3; o->s_gapo = 11; o->s_gape = 4;
	o->max_diff = 4; o->max_gapo = 1; o->max_gape = 6;
	o->report_best_only = 0;
	o->indel_end_skip = 5; o->max_del_occ = 10; o->max_entries = 2000000;
	o->n_batch = 0x40000;
	o->max_report_multi = 20;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	//o->seed_len = 32; o->max_seed_diff = 2;
	o->fnr = 0.06;
	o->n_threads = 1;
	o->max_top2 = 30;
	o->trim_qual = 0;
	o->rg = 0;
	o->verbose = 0;
	return o;
}

int bwa_cal_maxdiff(int l, double err, double thres)
{
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < 1000; ++k) {
		y *= l * err;
		x *= k;
		sum += elambda * y / x;
		if (1.0 - sum < thres) return k;
	}
	return 2;
}

// width must be filled as zero
static int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0; l = rbwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(rbwt, k - 1, l, c, &ok, &ol);
			k = rbwt->L2[c] + ok + 1;
			l = rbwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = rbwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

/*calculate the SA interval of one query sequence
 */
void jigsaw_cal_sa_reg_gap(bwt_t *const bwt[2], bwa_seq_t *seq, gap_stack_t *stack, const gap_opt_t *opt)
{
	bwt_width_t *w[2]; //, *seed_w[2];
	const ubyte_t *s[2];

	bwa_seq_t *p = seq;

	p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
	
	//sense_strand is init' here too: 0 +, 1 -, 2 .while opt->strand_mode 1 +, 2 -, 3 .
	p->sense_strand = 2; //opt->strand_mode - 1;
	
	s[0] = p->seq; s[1] = p->rseq;

	w[0] = w[1] = 0;
	w[0] = (bwt_width_t*)realloc(w[0], (p->len + 1) * sizeof(bwt_width_t));
	w[1] = (bwt_width_t*)realloc(w[1], (p->len + 1) * sizeof(bwt_width_t));
	memset(w[0], 0, (p->len + 1) * sizeof(bwt_width_t));
	memset(w[1], 0, (p->len + 1) * sizeof(bwt_width_t));


	//note that seq[0] is the forward sequence without reverse now, and seq[1] is the reverse complementary sequence
	//both need to be aligned to bwt[0], and bwt[1] should be used to calculate the bound

	bwt_cal_width(bwt[1], p->len, s[0], w[0]);
	bwt_cal_width(bwt[1], p->len, s[1], w[1]);


	// core function
	p->aln = bwt_match_gap(bwt[0], p->len, s, w, opt, &p->n_aln, stack);
	// store the alignment
	free(w[0]); free(w[1]);
}

bool jigsaw_check_seq_complexity (jigsaw_anchor_seq_t *anchor_seq)
{
	//if the seq is repetitive, return a 1, otherwise, return a 0
	bool is_repetitive = 1;
	ubyte_t *di_nt = (ubyte_t*) calloc (2, sizeof(ubyte_t) );
	di_nt[0] = anchor_seq->seq[0];
	di_nt[1] = anchor_seq->seq[1];
	for (int i = 2; i<anchor_seq->len; i++){
		if (anchor_seq->seq[i] != di_nt[i%2] ){
			is_repetitive = 0;
			break;
		}
	}
	free(di_nt);
	return is_repetitive;

}

/*
  this function returns a list of all hits of an anchor sequence
  for single anchor search
*/

void jigsaw_collect_anchor_hits (bwt_t *const bwt[2], jigsaw_anchor_seq_t *anchor_seq, const int *g_log_n,
                int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,
		int max_word_occ, const gap_opt_t *opt, list<jigsaw_word_hit_t*> *hits)
{
	//check if the anchor_seq is low-complexity, those seqs will highly affect the mapping speed. 
	//e.g. GTGTGTGTGT, ACACACACAC 
	if ( ( ! opt->allow_rep_anchor ) && jigsaw_check_seq_complexity (anchor_seq) ) return;
	gap_opt_t local_opt = *opt;
	
	//local_opt.max_diff = opt->max_word_diff;
	local_opt.max_diff = 0;//do not allow mismatch in single anchor search for now
	if( anchor_seq ->len >=18 && opt->max_diff >1) local_opt.max_diff = 1; // allow mismatch when the anchor is long
	//TODO: allow mismatch later, and set it as a parameter.  
	local_opt.max_gapo = local_opt.max_gape = 0;
	
	gap_stack_t *stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);
	// all the SA will be saved in anchor_seq->aln
	jigsaw_cal_sa_reg_gap (bwt, anchor_seq, stack, &local_opt);
	jigsaw_collect_word_hits_core (bwt[0], anchor_seq, hits);
	free(anchor_seq->aln);	
	gap_destroy_stack(stack);
}

/*
 * return a list of all hits of all non-repetitive words,
 * and also record the occurrence of each word in n_occ
 *
 * max_word_occ: the maximum number of occurrence allowed for a non-repetitive word
 * n_hits: return the number of hits found
 *
 */
void jigsaw_collect_word_hits (bwt_t *const bwt[2], jigsaw_word_t *words, int n_words, const int *g_log_n,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,
		int max_word_occ, const gap_opt_t *opt, list<jigsaw_word_hit_t*> *hits)
{
	int i;
	gap_opt_t local_opt = *opt;
	local_opt.max_diff = opt->max_word_diff;
	local_opt.max_gapo = local_opt.max_gape = 0;

	//the alignment of words need to be perfect, so we create a local (smaller) stack to
	//avoid unnecessarily big stack
	gap_stack_t *stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

	//find perfect match for each word

	//jigsaw_word_hit_t *hits = 0;
	int total_occ = 0;

	for (i = 0; i != n_words; ++i) {
		jigsaw_word_t *w = words+i;

		/*find all SA intervals stored at w->aln*/
		jigsaw_cal_sa_reg_gap (bwt, w, stack, &local_opt);

		/*the total number of hits of this word*/
		int word_occ = jigsaw_word_hits_total (w);
		w->n_occ = word_occ;

		if (word_occ <= max_word_occ) {
			//hits = (jigsaw_word_hit_t*)realloc (hits, (total_occ + word_occ) * sizeof (jigsaw_word_hit_t));

			/*calculate the query and target coordinates for all hits of the word, and save it to "hits"*/
			jigsaw_collect_word_hits_core (bwt[0], w, hits);

			total_occ += word_occ;
		}
	}

	//*n_hits = total_occ;
	gap_destroy_stack(stack);

	//return hits;
}


/*compare two hits by diagonal coordinates
 */

bool jigsaw_word_hits_comp_diagonal (const jigsaw_word_hit_t *a, const jigsaw_word_hit_t *b)
{
	int64_t da = a->pos_t - a->pos_q;
	int64_t db = b->pos_t - b->pos_q;
	return da < db ? true : false;
}


/* sort hits (ascendingly) by diagonal coordinates
 */
void jigsaw_sort_word_hits_by_diagonal (list<jigsaw_word_hit_t*> *hits)
{
	hits->sort (jigsaw_word_hits_comp_diagonal);
}


/* compare hits by strand, and then by diagonal
 */
int jigsaw_word_hits_comp_strand_diagonal (const void *a, const void *b)
{
	jigsaw_word_hit_t *ha = (jigsaw_word_hit_t*) a;
	jigsaw_word_hit_t *hb = (jigsaw_word_hit_t*) b;

	if (ha->strand > hb->strand) return 1;
	else if (ha->strand < hb->strand) return -1;
	else {//the same strand compare diagonal
		int64_t da = ha->pos_t - ha->pos_q;
		int64_t db = hb->pos_t - hb->pos_q;
		return da > db ? 1 : (da < db ? -1 : 0);
	}
}

/* sort hits (ascendingly) by diagonal coordinates
 */
void jigsaw_sort_word_hits_by_strand_diagonal (jigsaw_word_hit_t *hits, int n_hits)
{
	qsort (hits, n_hits, sizeof (jigsaw_word_hit_t), jigsaw_word_hits_comp_strand_diagonal);
}


/* compare hits by query coordinates
 */
bool jigsaw_word_hits_comp_pos_q (const jigsaw_word_hit_t *a, const jigsaw_word_hit_t *b)
{
	return a->pos_q < b->pos_q ? true : false;
}


/* sort hits (ascendingly) by query coordinates
 */
void jigsaw_sort_word_hits_by_pos_q (list<jigsaw_word_hit_t*> *hits)
{
	hits->sort (jigsaw_word_hits_comp_pos_q);
}



/* compare hits by strand, and then by pos_t
 */
bool jigsaw_word_hits_comp_strand_pos_t (const jigsaw_word_hit_t *a, const jigsaw_word_hit_t *b)
{
	//jigsaw_word_hit_t *ha = (jigsaw_word_hit_t*) a;
	//jigsaw_word_hit_t *hb = (jigsaw_word_hit_t*) b;

	if (a->strand < b->strand) return true;
	else if (a->strand > b->strand) return false;
	else {//the same strand compare pos_t
		return a->pos_t < b->pos_t ? true : false;
	}
}

/* sort hits (ascendingly) by diagonal coordinates
 */
void jigsaw_sort_word_hits_by_strand_pos_t (list<jigsaw_word_hit_t*> *hits)
{
	hits->sort(jigsaw_word_hits_comp_strand_pos_t);
	//qsort (hits, n_hits, sizeof (jigsaw_word_hit_t), jigsaw_word_hits_comp_strand_pos_t);
}


/*sort exons by target coordinates
 */
bool jigsaw_exon_comp_pos_t (const jigsaw_exon_t *a, const jigsaw_exon_t *b)
{
    //colinearity, descending
    if (a->colinear > b->colinear) return true;
    else if (a->colinear < b->colinear) return false;
    else return a->start_t < b->start_t ? true : false; //target coordinates, ascending
}


/* sort exons by target coordinates
 */

void jigsaw_sort_exons_by_pos_t (list<jigsaw_exon_t*> *exons)
{
    //sort all exons
    exons->sort (jigsaw_exon_comp_pos_t);
}



/*check if alignment of hits is colinear
 * hits must have been sorted according to pos_q (or equivalently wid)
 * so for a colinear clump, pos_t must be non-decreasing
 */
uint32_t jigsaw_hits_colinear (list<jigsaw_word_hit_t*> *hits)
{
	uint32_t colinear = 1;
	list<jigsaw_word_hit_t*>::iterator iter;
	iter = hits->begin();
	jigsaw_word_hit_t *a = *iter; ++iter;
	for (; iter != hits->end(); ++iter) {
		jigsaw_word_hit_t *b = *iter;
		if (a->pos_t > b->pos_t) {
			colinear = 0;
			break;
		}
		a = b;
	}
	return colinear;
}



/*calculate the coverage of hits in each exon
 * hits in the exon must have been sorted according to pos_q (or equivalently wid)
 */
float jigsaw_exon_coverage (jigsaw_exon_t *exon)
{
	int n_hits = exon->hits->size();
	jigsaw_word_hit_t *first, *last;
	list<jigsaw_word_hit_t *>::iterator iter = exon->hits->begin(); first = *iter;
	iter = exon->hits->end (); --iter; last = *iter;
	int n_words = last->wid - first->wid +1;
	return ((float)n_hits)/n_words;
}


/*1.for each exon, sort hits according to target coordinates
 *2. evaluate each exon in terms of colinearity, the number of hits, coverage, and uniqueness of words
 */
void jigsaw_eval_exons (list<jigsaw_exon_t*> *exons, int word_size)
{
	//int i;
	list<jigsaw_exon_t *>::iterator iter;

	for (iter = exons->begin(); iter != exons->end(); iter++) {
		jigsaw_exon_t *p = *iter;

		if (p->hits->size() == 0) continue; //this should never happen

		jigsaw_sort_word_hits_by_pos_q (p->hits);

		p->colinear = jigsaw_hits_colinear (p->hits);
		p->coverage = jigsaw_exon_coverage (p);
		//p->score = jigsaw_clump_uniqueness (p, words, l_pac); //log E-value of the clump

		jigsaw_word_hit_t *first = p->hits->front();
		jigsaw_word_hit_t *last = p->hits->back();

		p->strand = first->strand;
		p->start_q = first->pos_q;
		p->end_q = last->pos_q + word_size - 1;

		p->start_t = first->pos_t;
		p->end_t = last->pos_t + word_size - 1;

		//TODO: more properties?
	}
}


/* cluster hits according to diagonal coordinates
 * the idea borrowed from BLAT
 */
void jigsaw_group_hits_to_exons (list <jigsaw_word_hit_t*> *hits, int max_diff, int word_size, int n_words, list <jigsaw_exon_t *> *exons)
{
	int local_max_diff = int (max_diff);
        if( local_max_diff > 2 ) local_max_diff = 2;
	jigsaw_exon_t *curr_exon = 0;

	int64_t prev_hit_diag = 0, curr_hit_diag;
	uint32_t prev_strand = 2;
	//illegal value so that a new clump will be created at the very beginning

	list <jigsaw_word_hit_t *>::iterator iter;
	uint32_t current_exon_id = 0;

	//require each word to have unique occurrence in each cluster
	int* word_uniqueness = (int *)calloc(n_words, sizeof(int));
	for (iter = hits->begin(); iter != hits->end(); ++iter)
	{
		int wid = abs ((*iter)->wid);
		word_uniqueness[wid]++;
	}


	//fprintf (stderr, "group exons ...\n");
	for (iter = hits->begin(); iter != hits->end(); ++iter)
	{
		jigsaw_word_hit_t *p = *iter;
		if (word_uniqueness[abs(p->wid)]>1)
			continue;

		curr_hit_diag = p->pos_t - p->pos_q;
		if (p->strand != prev_strand || curr_hit_diag - prev_hit_diag > local_max_diff) {
			//create a new exon

			//initialization
			curr_exon = (jigsaw_exon_t *) calloc (1, sizeof (jigsaw_exon_t));
			curr_exon->n_mm = curr_exon->n_gapo_t = curr_exon->n_gapo_q = curr_exon->n_gape_t = curr_exon->n_gape_q = 0;

			curr_exon->exon_id = current_exon_id++;
			//fprintf(stderr, "create a new exon id=%d\n", curr_exon->exon_id);

			curr_exon->hits = new list<jigsaw_word_hit_t*>;
			curr_exon->is_first = curr_exon->is_last = 1;
			exons->push_back (curr_exon);
		}

		//add the hits to the current clump
		//double the capacity of the array if necessary
		curr_exon->hits->push_back (p);
		//fprintf (stderr, "exon_id=%d, start_t=%d, start_q=%d\n", curr_exon->exon_id, p->pos_t, p->pos_q);

		prev_hit_diag = curr_hit_diag;
		prev_strand = p->strand;
	}
	free (word_uniqueness);

	//fprintf (stderr, "%d exons identified\n", exons->size());

	//evaluate exons: calculate colinearity, etc
	jigsaw_eval_exons (exons, word_size);
}


/*copy one clump
 */
/*
void jigsaw_copy_clump (jigsaw_clump_t *to, const jigsaw_clump_t *from)
{
	*to = *from;
	to->hits = calloc (from->n_hits, sizeof (jigsaw_word_hit_t));
	memcpy (to->hits, from->hits, from->n_hits * sizeof (jigsaw_word_hit_t));
}
*/

/*filter clumps according to colinearity and score
 * colinear clumps with a score <= max_score is copied to a new array and returned
 */

/*
jigsaw_clump_t* jigsaw_filter_clumps (jigsaw_clump_t* clumps, int n_clumps, float max_score, int *n_filtered_clumps)
{
	int i, max_clumps = 4;
	jigsaw_clump_t *filtered_clumps = (jigsaw_clump_t*) calloc (max_clumps, sizeof (jigsaw_clump_t));

	int _n_filtered_clumps = 0;

	for (i = 0; i != n_clumps; ++i)
	{
		jigsaw_clump_t *p = clumps + i;
		if (p->colinear != 1 || p->score > max_score) continue;

		if (max_clumps <= _n_filtered_clumps) { //double the capacity of the array if necessary
			max_clumps <<= 1;
			filtered_clumps = (jigsaw_clump_t*) realloc (filtered_clumps, max_clumps * sizeof (jigsaw_clump_t));
		}

		jigsaw_copy_clump (filtered_clumps+_n_filtered_clumps, p);
		++_n_filtered_clumps;
	}

	*n_filtered_clumps = _n_filtered_clumps;

	return filtered_clumps;
}
*/



/*extend a hit on both ends by exact matches
 * direction: 0-extend on the right, 1-extend on the left, 2-extend both side
 */

void jigsaw_extend_exon_exact (jigsaw_exon_t *exon, bwa_seq_t *seq, int direction,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac)
{
	ubyte_t *query_seq = exon->strand ? seq->rseq : seq->seq;

	jigsaw_exon_t *p = exon;


	//extend on the left
	if (direction != 0) {
		int i = 1;
		while (p->start_q -i >= 0 && p->start_t -i >= 0) {
			int64_t kq = p->start_q -i, kt = p->start_t -i;

			ubyte_t q = query_seq [kq];
			ubyte_t t = get_pacseq_base (pacseq, kt);
			//ubyte_t t = pacseq[kt>>2] >> ((~kt&3)<<1) & 3;
			if (q != t) break;

			++i;
		}
		--i;
		p->start_q -= i; p->start_t -= i;
	}

	//extend on the right
	if (direction != 1) {
		int i = 1;
		while (p->end_q +i < seq->len && p->end_t +i < l_pac) {
			int64_t kq = p->end_q +i, kt = p->end_t +i;

			ubyte_t q = query_seq [kq];
			ubyte_t t = get_pacseq_base (pacseq, kt);
			//ubyte_t t = pacseq[kt>>2] >> ((~kt&3)<<1) & 3;
			if (q != t) break;

			++i;
		}
		--i;
		p->end_q += i; p->end_t += i;
	}
}


/*extend a hit on both ends by inexact matches allowing for substitutions but not indels
 * direction: 0-extend on the right, 1-extend on the left, 2-extend both side
 */

void jigsaw_extend_exon_inexact (jigsaw_exon_t *exon, bwa_seq_t *seq, int direction,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac, bool stop_if_neg_score, int opt_max_diff)
{
	ubyte_t *query_seq = exon->strand ? seq->rseq : seq->seq;

	jigsaw_exon_t *p = exon;
	int match = 1, mismatch = -3, max_diff = 2;
	if (opt_max_diff <2) max_diff = opt_max_diff;
	//TODO: implicitly require 4 matches at the end

	//extend on the left
	if (direction != 0) {
		int i = 0, end = 0;
		int score = 0, max_score = 0, diff = 0, min_diff = max_diff;

		while (p->start_q -i >= 0 && p->start_t -i >= 0) {
			int64_t kq = p->start_q -i, kt = p->start_t -i;

			ubyte_t q = query_seq [kq];
			ubyte_t t = get_pacseq_base (pacseq, kt);
			//ubyte_t t = pacseq[kt>>2] >> ((~kt&3)<<1) & 3;
			if (q == t) score += match;
			else {
				score += mismatch;
				diff++;
				if (diff > max_diff) break;
			}

			if (stop_if_neg_score && score < 0) break;

			//if the score is not the optimal, but it can extend to the very end of the read
			//we still want to accept the extension
			if (score > max_score || i == p->start_q) {
				max_score = score; end = i; min_diff = diff;
			}
			++i;
		}
		//if (end > 0) --end; //extendible
		p->start_q -= end; p->start_t -= end;
		p->n_mm += min_diff;
	}

	//extend on the right
	if (direction != 1) {
		int i = 0, end = 0;
		int score = 0, max_score = 0, diff = 0, min_diff = max_diff;

		while (p->end_q +i < seq->len && p->end_t +i < l_pac) {
			int64_t kq = p->end_q +i, kt = p->end_t +i;

			ubyte_t q = query_seq [kq];
			//ubyte_t t = pacseq[kt>>2] >> ((~kt&3)<<1) & 3;
			ubyte_t t = get_pacseq_base (pacseq, kt);
			if (q == t) score += match;
			else {
				score += mismatch;
				diff++;
				if (diff > max_diff) break;
			}

			if (stop_if_neg_score && score < 0) break;

			//if the score is not the optimal, but it can extend to the very end of the read
			//we still want to accept the extension
			if (score > max_score || p->end_q+i == seq->len - 1) {
				max_score = score; end = i; min_diff = diff;
			}
			++i;
		}
		//if (end > 0) --end; //extendible
		p->end_q += end; p->end_t += end;
		p->n_mm += min_diff;
	}
}


void jigsaw_extend_exons (list<jigsaw_exon_t*> *exons, bwa_seq_t *seq,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac, uint32_t exact, bool stop_if_neg_score, int opt_max_diff)
{
	list<jigsaw_exon_t*>::iterator iter;
	for (iter = exons->begin(); iter != exons->end (); ++iter) {
		jigsaw_exon_t *p = *iter;
		exact ? jigsaw_extend_exon_exact (p, seq, 2, l_pac, pacseq, ntpac) : jigsaw_extend_exon_inexact (p, seq, 2, l_pac, pacseq, ntpac, stop_if_neg_score, opt_max_diff);
	}
}




/* fill in gaps in an exon to the reference genome
 * TODO: now the mismataches are counted only in the gaps, but now in the words themselves
 */

void jigsaw_refine_exon_aln_core (jigsaw_exon_t *exon, bwa_seq_t *seq,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac, int band_width)
{
	int path_len;
	jigsaw_exon_t *p = exon;

	AlnParam ap = aln_param_bwa;
	//TODO: need more sophisticated rules; the bound is too loose here
	ap.band_width = band_width>1 ? band_width: 1;
	//TODO, band_width > 1 for now, avoid crashes in aln_global_core
	
	//if (strcmp (seq->name, "HWI-ST155_05162111902026#0") == 0)
	//{
	//	fprintf (stderr, "found\n");
	//}


	//initialization
	//p->n_mm = p->n_gapo_t = p->n_gapo_q = p->n_gape_t = p->n_gape_q = 0;

	//the line above is unnessary now since there is inexact extending step -j
	//need to find a way to sum up the mismatches...
	if (p->hits->size() < 2) return;

	//fprintf (stderr, "refine exons ...\n");

	list<jigsaw_word_hit_t *>::iterator prev, next;
	prev = next = p->hits->begin(); ++next;
	for (; next != p->hits->end (); ++prev, ++next) {

		int64_t ref_start = (*prev)->pos_t + seq->word_size;
		int ref_len = (*next)->pos_t - ref_start;

		int start = (*prev)->pos_q + seq->word_size;
		int len = (*next)->pos_q - start;

		//fprintf (stderr, "exon_id=%d, start_t1=%d, start_q1=%d, start_t2=%d, start_q2=%d\n", exon->exon_id, (*prev)->pos_t, (*prev)->pos_q, (*next)->pos_t, (*next)->pos_q);

		//TODO: double check if there are bugs in the following lines in rare cases
		if (ref_len <= 0 && len <= 0 && len == ref_len ) continue;
		else if (ref_len <= 0 && len > ref_len ) {
			//no nucleotide in reference
			p->n_gapo_t += 1; p->n_gape_t += (len - ref_len - 1);
			continue;
		}
		else if (len <= 0 && ref_len > len) {
			//no nucleotide in query
			p->n_gapo_q += 1; p->n_gape_q += (ref_len - len - 1);
			continue;
		}

		//both are not empty
		ubyte_t* ref_seq = (ubyte_t*)calloc(ref_len, sizeof(ubyte_t));

		int64_t k, l;
		for (k = ref_start, l= 0; k < ref_start + ref_len; ++k)
			ref_seq[l++] = get_pacseq_base (pacseq, k); //pacseq[k>>2] >> ((~k&3)<<1) & 3;

		ubyte_t *query_seq = (exon->strand ? seq->rseq : seq->seq) + start;


		path_t *path = (path_t*)calloc(ref_len+len, sizeof(path_t));
		aln_global_core(ref_seq, ref_len, query_seq, len, &ap, path, &path_len);

		//calculate the number of mismatches and gaps
		uint32_t n_mm, n_gapo_t, n_gapo_q, n_gape_t, n_gape_q;
		jigsaw_cal_diff (path, path_len, ref_seq, ref_len, query_seq, len,
			&n_mm, &n_gapo_t, &n_gapo_q, &n_gape_t, &n_gape_q);

		free (ref_seq); free (path);

		p->n_mm += n_mm;
		p->n_gapo_t += n_gapo_t; p->n_gapo_q += n_gapo_q;
		p->n_gape_t += n_gape_t; p->n_gape_q += n_gape_q;
	}
}


void jigsaw_refine_exon_aln (list<jigsaw_exon_t*> *exons, bwa_seq_t *seq,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac, int band_width)
{
	list<jigsaw_exon_t*>::iterator iter;
	for (iter = exons->begin(); iter != exons->end(); ++iter) {
		jigsaw_refine_exon_aln_core (*iter, seq, l_pac, pacseq, ntpac, band_width);
	}
}



/*an approximate estimate of E-value
 * count each word only once if multiple hits exists in the alignment
 */
float jigsaw_spliced_aln_cluster_uniqueness (jigsaw_spliced_aln_cluster_t *aln,
		jigsaw_word_t *words, int n_words, int64_t l_pac)
{
	float log_occ = 0;
	uint32_t *u = (uint32_t *) calloc (n_words, sizeof(uint32_t));
	int n_hits_uniq = 0;

	list <jigsaw_word_hit_t *>::iterator iter;
	for (iter = aln->hits->begin(); iter != aln->hits->end(); ++iter) {
		jigsaw_word_hit_t *h = *iter;
		int wid = abs(h->wid);
		if (u[wid]) continue;

		++u[wid];
		++n_hits_uniq;
		jigsaw_word_t* w = words + wid;
		log_occ += log10((float)w->n_occ);
	}
	free (u);
	return log_occ - (n_hits_uniq -1) * log10 ((float)l_pac);
}


/*evaluate candidate alignment according to uniqueness of words
 */

static void jigsaw_eval_spliced_aln_clusters (list<jigsaw_spliced_aln_cluster_t *> *clusters,
		jigsaw_word_t *words, int n_words, int word_size, int64_t l_pac)
{
	list<jigsaw_spliced_aln_cluster_t *>::iterator iter;

	for (iter = clusters->begin(); iter != clusters->end(); ++iter) {
		jigsaw_spliced_aln_cluster_t *p = *iter;

		p->score = jigsaw_spliced_aln_cluster_uniqueness (p, words, n_words, l_pac); //log E-value of the clump

		//p->strand = p->hits[0].strand;
		//p->start_q = p->hits[0].pos_q;
		//p->end_q = p->hits[p->n_hits-1].pos_q + word_size - 1;

		//p->start_t = p->hits[0].pos_t;
		//p->end_t = p->hits[p->n_hits-1].pos_t + word_size - 1;

		//TODO: more properties?
	}
}

/*
 * not in use
 
uint32_t jigsaw_spliced_aln_cluster_unique_hits (jigsaw_spliced_aln_cluster_t *aln, int n_words)
{

	uint32_t *u = (uint32_t *) calloc (n_words, sizeof(uint32_t));
	uint32_t n_hits_uniq = 0;

	list <jigsaw_word_hit_t *>::iterator iter;
	for (iter = aln->hits->begin(); iter != aln->hits->end(); ++iter) {
		jigsaw_word_hit_t *h = *iter;
		int wid = abs(h->wid);
		if (u[wid]) continue;

		++u[wid];
		++n_hits_uniq;
	}
	free (u);
	return n_hits_uniq;
}
*/

/*hits should have been sorted by strand and then by pos_t before calling this function*/
void jigsaw_group_hits_to_spliced_aln_clusters (list<jigsaw_word_hit_t*> *hits, int max_intron_size,
		jigsaw_word_t *words, int n_words, int word_size, int64_t l_pac, list<jigsaw_spliced_aln_cluster_t *> *clusters)
{
	jigsaw_spliced_aln_cluster_t *curr_clust = 0;

	int64_t prev_pos_t = 0, curr_pos_t;
	uint32_t prev_strand = 2;
	//illegal value so that a new clump will be created at the very beginning
	int prev_wid_direction = 0, curr_wid_direction = 0; //this is to make sure hits->wid in the cluster are in the same direction, 0 means not determined.
	int prev_wid = -1, curr_wid;

	list<jigsaw_word_hit_t *>::iterator iter;
	uint32_t curr_cluster_id = 0;
	for (iter = hits->begin (); iter != hits->end(); ++iter)
	{
		jigsaw_word_hit_t *p = *iter;
		curr_pos_t = p->pos_t;
		curr_wid = abs(p->wid);
		if(prev_wid != -1)
		{
			int wid_diff = curr_wid - prev_wid;
			if ( wid_diff != 0) {curr_wid_direction = abs( wid_diff)/ wid_diff;}
			else {curr_wid_direction = 2;} // 2 means the same direction
		}


		if (p->strand != prev_strand || curr_pos_t - prev_pos_t > 2 * max_intron_size || 
			curr_wid_direction ==2 	|| 
			(prev_wid_direction != 0 && prev_wid_direction != curr_wid_direction ) ) {
		    //2 times maxintron size because of possible inner exons
		    //TODO: after increasing this, we missed some alignments, maybe because the uniqueness issue

/*
			if (curr_clust && jigsaw_spliced_aln_cluster_unique_hits (curr_clust, n_words) < (uint32_t) n_words/2)
			{
				if (curr_clust->hits) delete hits;
				//free (curr_clust);
				clusters->pop_back();
			}
*/
			//prev_wid_direction = 0;
			curr_wid_direction = 0; 
			//prev_wid = -1;
			//curr_wid = -1;

			//create a new candidate alignment

			curr_clust = (jigsaw_spliced_aln_cluster_t*) calloc (1, sizeof (jigsaw_spliced_aln_cluster_t));
			curr_clust->cluster_id = curr_cluster_id++;

			//initialization
			//curr_clust->n_junctions = 0;
			curr_clust->junctions = 0;
			curr_clust->hits = new list<jigsaw_word_hit_t *>;
			//add to the list
			clusters->push_back (curr_clust);
		}

		//add the hit to the current alignment
		//double the capacity of the array if necessary

		curr_clust->hits->push_back (p);
		prev_wid_direction = curr_wid_direction;

		prev_pos_t = p->pos_t;
		prev_wid = curr_wid;
		prev_strand = p->strand;
	}

	//*n_clusters = _n_clusters;

	//evaluate aln: calculate start, end and score
	jigsaw_eval_spliced_aln_clusters (clusters, words, n_words, word_size, l_pac);
}


/*compare spliced alignment by score (ascending)
 */

bool jigsaw_spliced_aln_cluster_comp_score (const jigsaw_spliced_aln_cluster_t *a, const jigsaw_spliced_aln_cluster_t *b)
{
	return a->score < b->score ? true : false;
}


/*compare spliced alignment by score (ascending)
 *
 */
void jigsaw_sort_spliced_aln_cluster_by_score (list<jigsaw_spliced_aln_cluster_t*> *clusters)
{
	//sort all clumps
	clusters->sort (jigsaw_spliced_aln_cluster_comp_score);
}



/*release the memory of spliced alignment
 */

void jigsaw_destroy_word_hits (list<jigsaw_word_hit_t*> *hits)
{
	list<jigsaw_word_hit_t*>::iterator iter;
	for (iter = hits->begin(); iter != hits->end(); ++iter) {
		jigsaw_word_hit_t *p = *iter;
		free(p);
		//no need to remove each hit here
	}
	//delete (hits);
}


void jigsaw_destroy_exons (list<jigsaw_exon_t*> *exons)
{
	list<jigsaw_exon_t*>::iterator iter;
	for (iter = exons->begin(); iter != exons->end(); ++iter) {
		jigsaw_exon_t *p = *iter;
		delete (p->hits); free(p);
		//no need to remove each hit here
	}
	//delete (exons);
}

void jigsaw_destroy_junctions (list<jigsaw_junction_t*> *junctions)
{
	list<jigsaw_junction_t*>::iterator iter;
	for (iter = junctions->begin(); iter != junctions->end(); ++iter) {
		jigsaw_junction_t *p = *iter;
		free(p);
		//no need to remove each hit here
	}
	//delete (junctions);
}

void jigsaw_destroy_spliced_aln_clusters (list<jigsaw_spliced_aln_cluster_t*> *clusters)
{
	list<jigsaw_spliced_aln_cluster_t*>::iterator iter;
	for (iter = clusters->begin(); iter != clusters->end(); ++iter) {
		jigsaw_spliced_aln_cluster_t *p = *iter;

		if (p->hits) delete (p->hits); //no need to delete each hit here
		if (p->exons) {
			jigsaw_destroy_exons (p->exons); delete (p->exons);
		}
		if (p->junctions) {
			jigsaw_destroy_junctions (p->junctions); delete (p->junctions);
		}
		free (p);
	}
	//delete (clusters);
}

void jigsaw_destroy_spliced_aln (list<jigsaw_spliced_aln_t*> *aln)
{
	list<jigsaw_spliced_aln_t*>::iterator iter;
	for (iter = aln->begin(); iter != aln->end(); ++iter) {
		jigsaw_spliced_aln_t *p = *iter;
		if(p->junctions != NULL) delete (p->junctions);
		free (p);
	}
	//delete (aln);
}

void jigsaw_get_partner_splice_site (uint32_t type, bool sense_strand, ubyte_t *ss)
//here the type is the original site, so if it is Donor, the partner is accepter
{
//	ubyte_t *ss = (ubyte_t*) calloc (2, sizeof(ubyte_t) );
	if(type == SPLICE_DONOR ) {
		if (sense_strand) {
			ss[0] = BASE_C; ss[1] = BASE_T;
		}
		else {
			ss[0] = BASE_A; ss[1] = BASE_G;
		}
	}
	else {
		if (sense_strand) {
			ss[0] = BASE_A; ss[1] = BASE_C;
		}
		else {
			 ss[0] = BASE_G; ss[1] = BASE_T;
		}
	}
//	return ss;
}


inline bool jigsaw_is_donor_splice_site (ubyte_t *ss, bool revcom)
{
	if (revcom)
		return ss[0] == BASE_A && ss[1] == BASE_C ? true : false;
	else
		return ss[0] == BASE_G && ss[1] == BASE_T ? true : false;
}

inline bool jigsaw_is_acceptor_splice_site (ubyte_t *ss, bool revcom)
{
	if (revcom)
		return ss[0] == BASE_C && ss[1] == BASE_T ? true : false;
	else
		return ss[0] == BASE_A && ss[1] == BASE_G ? true : false;
}


inline bool jigsaw_is_splice_site (ubyte_t *ss, uint32_t type, bool revcom)
{
	if (type == SPLICE_DONOR)
		return jigsaw_is_donor_splice_site (ss, revcom);
	else if (type == SPLICE_ACCEPTOR)
		return jigsaw_is_acceptor_splice_site (ss, revcom);
	else
		return false;
}


//indicate if the current position (the first or last base of an intron) is a splice site
inline bool jigsaw_is_splice_site_pos (uint64_t pos, uint32_t type, bool sense_strand, int64_t l_pac, const ubyte_t *pacseq)
{
	ubyte_t ss[2];
	if ((sense_strand == 0 && type == SPLICE_DONOR)||
		(sense_strand == 1 && type == SPLICE_ACCEPTOR))
	{
		ss[0] = get_pacseq_base (pacseq, pos);
		ss[1] = get_pacseq_base (pacseq, pos+1);
	}
	else
	{
		ss[0] = get_pacseq_base (pacseq, pos-1);
		ss[1] = get_pacseq_base (pacseq, pos);
	}
	return jigsaw_is_splice_site (ss, type, sense_strand);
}

/*
 * not in use anymore
 * extract splice site sequences at the position
  * direction: 0, downstream; 1, get upstream
 * sense_strand: 0 for "+" 1 for "-"
 */
/*
void jigsaw_enum_partner_splice_site_denovo (int64_t pos, int direction, int max_intron_size, int min_intron_size,
		uint32_t sense_strand, int64_t l_pac, const ubyte_t *pacseq, list<splice_site_t*>* splice_sites)
{
	ubyte_t ss[2];
	//int min_intron_size = 20;

	if (direction == 0)
	{
		//look for downstream introns
		ss[0] = get_pacseq_base (pacseq, pos);
		ss[1] = get_pacseq_base (pacseq, pos+1);

		//look for a donor if this is + strand, or acceptor if this is - strand
		if (sense_strand == 0)
		{
			if (jigsaw_is_donor_splice_site (ss, 0) == false) return;
		}
		else
		{
			if (jigsaw_is_acceptor_splice_site (ss, 1) == false) return;
		}

		for (int64_t iend = pos + min_intron_size - 1; iend <= pos + max_intron_size - 1 && iend < l_pac; iend++)
		{
			ss[0] = get_pacseq_base (pacseq, iend-1);
			ss[1] = get_pacseq_base (pacseq, iend);

			//look for an acceptor if this is + strand, or donor if this is - strand
			if (sense_strand == 0)
			{
				if (jigsaw_is_acceptor_splice_site (ss, 0) == false) continue;
			}
			else
			{
				if (jigsaw_is_donor_splice_site (ss, 1) == false) continue;
			}

			splice_site_t *p = (splice_site_t *) calloc(1, sizeof (splice_site_t));
			p->type = sense_strand == 0 ? SPLICE_ACCEPTOR : SPLICE_DONOR;
			p->strand = sense_strand;
			p->pos = iend;
			p->n_partners = p->max_partners = 0;
			p->partners = 0;
			splice_sites->push_back (p);
		}
	}
	else
	{
		//look for upstream introns
		ss[0] = get_pacseq_base (pacseq, pos-1);
		ss[1] = get_pacseq_base (pacseq, pos);

		//look for an acceptor if this is + strand, or donor if this is - strand
		if (sense_strand == 0)
		{
			if (jigsaw_is_acceptor_splice_site (ss, 0) == false) return;
		}
		else
		{
			if (jigsaw_is_donor_splice_site (ss, 1) == false) return;
		}

		for (int64_t istart = pos - min_intron_size + 1; istart >= pos - max_intron_size + 1 && istart >= 0; istart--)
		{
			ss[0] = get_pacseq_base (pacseq, istart);
			ss[1] = get_pacseq_base (pacseq, istart+1);

			//look for an donor if this is + strand, or acceptor if this is - strand
			if (sense_strand == 0)
			{
				if (jigsaw_is_donor_splice_site (ss, 0) == false) continue;
			}
			else
			{
				if (jigsaw_is_acceptor_splice_site (ss, 1) == false) continue;
			}

			splice_site_t *p = (splice_site_t *) calloc(1, sizeof (splice_site_t));
			p->type = sense_strand == 0 ? SPLICE_DONOR : SPLICE_ACCEPTOR;
			p->strand = sense_strand;
			p->pos = istart;
			p->n_partners = p->max_partners = 0;
			p->partners = 0;
			splice_sites->push_back (p);
		}
	}
}
*/


/*
 * this function is not used anymore
 * extract splice sites in the interval start and end
 * sense_strand: 0 for "+" 1 for "-"
 * type: SPLICE_DONOR, SPLICE_ACCEPTOR
 * direction: from right to left (1) or from left to right (0)
 */
 /*
void jigsaw_enum_splice_site_denovo (int64_t start, int64_t end, uint32_t sense_strand, uint32_t type, int direction,
		int64_t l_pac, const ubyte_t *pacseq, list<splice_site_t*>* splice_sites)
{
	xassert (type == SPLICE_DONOR || type == SPLICE_ACCEPTOR, "incorrect splice site");

	ubyte_t ss[2];
	//int min_intron_size = 20;
	for (int64_t pos0 = start; pos0 < end-1; ++pos0)
	{
		int64_t pos = direction ? end -2 + start - pos0 : pos0;

		//look for downstream introns
		ss[0] = get_pacseq_base (pacseq, pos);
		ss[1] = get_pacseq_base (pacseq, pos+1);

		if (jigsaw_is_splice_site (ss, type, sense_strand))
		{
			splice_site_t *p = (splice_site_t *) calloc(1, sizeof (splice_site_t));
			p->type = type;
			p->strand = sense_strand;
			p->pos = pos;

			if ((sense_strand == 0 && type == SPLICE_ACCEPTOR)||
			   (sense_strand == 1 && type == SPLICE_DONOR)) p->pos++; //always point to the base next to the exon

			p->n_partners = p->max_partners = 0;
			p->partners = 0;
			splice_sites->push_back (p);
		}
	}
}
*/

int jigsaw_locate_junc_one_anchor_with_anno_downstream(jigsaw_exon_t *exon, bwa_seq_t *seq, int64_t min_ue_end_t, int64_t max_ue_end_t, int n_backsearch,int64_t l_pac,  const ubyte_t *pacseq, const gap_opt_t *opt, list<jigsaw_junction_t*> *junctions, list<jigsaw_exon_t*> *exons)
{
    int local_max_diff = (int) opt->max_diff;
    if(local_max_diff >2 ) local_max_diff = 2;
    //int best_diff = opt->max_diff+1;
    int best_diff = local_max_diff;
    
    int num_junc_found_in_anno = 0;
    int64_t k;
    jigsaw_junction_t *junction = 0;
    jigsaw_exon_t *de = 0;
    //analog to the two_anchor search
    AlnParam ap = aln_param_bwa;
    //ap.band_width = opt->max_diff;
    ap.band_width = local_max_diff >1 ? local_max_diff: 1;
  
    int gap_start_q = exon->end_q - n_backsearch +1;
    int gap_len = seq->len - gap_start_q ;
  
    ubyte_t *query_seq = exon->strand ? seq->rseq : seq->seq;
    
    ubyte_t *gap_seq = query_seq + gap_start_q;
    ubyte_t *ref_seq = (ubyte_t *) malloc(gap_len * sizeof(ubyte_t));
    
    uint32_t *adjust_diff =new uint32_t[ n_backsearch +1 ];
    int diff_track = 0;
    for ( int ki =0; ki< n_backsearch + 1; ki++) {
	    int64_t kt =  exon->end_t - ki;
	    int64_t kq = exon->end_q - ki;
	    ubyte_t t = get_pacseq_base (pacseq, kt);
	    ubyte_t q = query_seq [kq];
	    if (q != t) diff_track++;
	    adjust_diff[ki] = diff_track;
     }
    //sense_strand: 0 is "+", 1 is "-"
    for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
	    
	    if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != exon->strand ) ) continue;
	    if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == exon->strand ) ) continue;
		    
	  int64_t ue_end_q = exon->end_q - n_backsearch ;
	  uint32_t type = sense_strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;
	  for (int64_t ue_end_t = min_ue_end_t; ue_end_t <= max_ue_end_t; ++ue_end_t, ue_end_q++) {
	      if (seq->len - ue_end_q -1 < opt->known_junc_min_anchor ) break;
	      int64_t intron_start_t = ue_end_t + 1;
	      splice_site_t* known_ss= retrieve_splice_site (opt->splice_site_map, intron_start_t, sense_strand, type);
	      if(!known_ss) continue;
	      //check all the partner sites recorded in splice_site_map
	      int i=0;
	      //put the first half of the sequence into ref_seq
	      for (k = min_ue_end_t + 1; k <= ue_end_t; ++k, ++i)
		      ref_seq[i] = get_pacseq_base (pacseq,k);
		      
	      for(int partner_i = 0; partner_i < known_ss->n_partners; partner_i++ ){
		      splice_site_t* partner_ss = known_ss->partners[partner_i];
		      int64_t de_start_t = partner_ss->pos +1;
		      if (de_start_t + seq->len - ue_end_q - 2 >= l_pac) continue;
		      //avoid exceed reference boundary
		      //start to do alignment near this site
		      int ii=i;
		      //put the other half into ref_seq
		      for (k = de_start_t; ii<gap_len; ++ii, ++k)
			      ref_seq[ii] = get_pacseq_base (pacseq,k);
		      path_t *path = (path_t*)calloc(gap_len * 2, sizeof(path_t));
		      int path_len;
		      aln_global_core(ref_seq, gap_len, gap_seq, gap_len, &ap, path, &path_len);
		      //calculate the number of mismatches and gaps
		      uint32_t n_mm, n_gapo_t, n_gapo_q, n_gape_t, n_gape_q, an_mm;
		      jigsaw_cal_diff (path, path_len, ref_seq, gap_len, gap_seq, gap_len,
		      &n_mm, &n_gapo_t, &n_gapo_q, &n_gape_t, &n_gape_q);
		      free (path);
		      an_mm =0;
		      if (ue_end_t - exon->end_t <=0)  an_mm =  adjust_diff[ exon->end_t-ue_end_t ];
		      int diff = n_mm + n_gapo_t + n_gapo_q + n_gape_t + n_gape_q - an_mm;
		      //better junction found
		      if(opt->report_best_only){
			  if (diff <= local_max_diff && diff < best_diff) {
				   
				   best_diff = diff;
				   
				  //create a de first
				  if (!de) de = (jigsaw_exon_t*)calloc(1, sizeof(jigsaw_exon_t));
				  de->n_mm= de->n_gapo_t = de->n_gape_t = de->n_gapo_q = de->n_gape_q = 0;
				  de->colinear = 1; de->coverage = 1;
				  de->strand = exon->strand;
				  de->start_t = de_start_t;
				  de->end_t = de_start_t + seq->len - ue_end_q - 2;
				  de->start_q = ue_end_q + 1;
				  de->end_q = seq->len -1;
				  de->hits = 0;
				  de->is_first = 0;
				  de->is_last = 1;

				  //modify the upstream exon
				  //exon->end_q = ue_end_q;
				  //exon->end_t = ue_end_t;
				  exon->is_last = 0;
				  
				  //create a junction
				  if (!junction) junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
				  junction->uexon = exon;
				  junction->dexon = de;
				  junction->sense_strand = sense_strand;
				  junction->n_mm = n_mm;
				  junction->an_mm = an_mm;
				  junction->n_gapo_t = n_gapo_t;
				  junction->n_gapo_q = n_gapo_q;
				  junction->n_gape_t = n_gape_t;
				  junction->n_gape_q = n_gape_q;
				  junction->start_q = ue_end_q;
				  junction->end_q = junction->start_q + 1;
				  junction->start_t = ue_end_t;
				  junction->end_t = de->start_t;
				  junction->logistic_prob = 1;
				   // junction->logistic_prob = get_splice_score(pacseq, ue_end_t, de->start_t);
			    }
		      }
		      else{
			 //push all the results into the lists
			  if (diff <= local_max_diff ) {
			  
			     //create a de first
			     de = (jigsaw_exon_t*)calloc(1, sizeof(jigsaw_exon_t));
			     de->n_mm= de->n_gapo_t = de->n_gape_t = de->n_gapo_q = de->n_gape_q = 0;
			     de->colinear = 1; de->coverage = 1;
			     de->strand = exon->strand;
			     de->start_t = de_start_t;
			     de->end_t = de_start_t + seq->len - ue_end_q - 2;
			     de->start_q = ue_end_q + 1;
			     de->end_q = seq->len -1;
			     de->hits = 0;
			     de->is_first = 0;
			     de->is_last = 1;

			     //modify the upstream exon
			     //exon->end_q = ue_end_q;
			     //exon->end_t = ue_end_t;
			     exon->is_last = 0;
			     
			     //create a junction
			     junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
			     junction->uexon = exon;
			     junction->dexon = de;
			     junction->sense_strand = sense_strand;
			     junction->n_mm = n_mm;
			     junction->an_mm = an_mm;
			     junction->n_gapo_t = n_gapo_t;
			     junction->n_gapo_q = n_gapo_q;
			     junction->n_gape_t = n_gape_t;
			     junction->n_gape_q = n_gape_q;
			     junction->start_q = ue_end_q;
			     junction->end_q = junction->start_q + 1;
			     junction->start_t = ue_end_t;
			     junction->end_t = de->start_t;
			     //junction->logistic_prob = get_splice_score(pacseq, ue_end_t, de->start_t);
			     junction->logistic_prob = 1;
			     junctions->push_back(junction);
			     exons->push_back(de);
			     num_junc_found_in_anno++;
			     				
			    }
		      }

	      }

	  }
	  
  }
		
  free(ref_seq);
  if(opt->report_best_only && junction){
      junctions->push_back(junction);
      exons->push_back(de);
      num_junc_found_in_anno++;
  }
  delete[] adjust_diff;
  return num_junc_found_in_anno;
}


int jigsaw_locate_junc_one_anchor_denovo_downstream(jigsaw_exon_t *exon, bwa_seq_t *seq, int64_t min_ue_end_t, int64_t max_ue_end_t, int n_backsearch, int64_t l_pac,const ubyte_t *pacseq, const ubyte_t *ntpac, bwt_t *const bwt[2], const int *g_log_n, const gap_opt_t *opt, list<jigsaw_junction_t*> *junctions, list<jigsaw_exon_t*> *exons)
{
	if (!  (seq->len - exon->end_q - 1 + opt->max_overhang >= opt->min_anchor && seq->len - exon->end_q - 1 < 2* seq->word_size) ) return 0;
    int num_junc_found_denovo = 0;
    ubyte_t *query_seq = exon->strand ? seq->rseq : seq->seq;
    jigsaw_exon_t *de = 0;
    jigsaw_junction_t *junction = 0;
    
    uint32_t *adjust_diff =new uint32_t[ n_backsearch +1 ];
    int diff_track = 0;
    for ( int ki =0; ki< n_backsearch + 1; ki++) {
	    int64_t kt =  exon->end_t - ki;
	    int64_t kq = exon->end_q - ki;
	    ubyte_t t = get_pacseq_base (pacseq, kt);
	    ubyte_t q = query_seq [kq];
	    if (q != t) diff_track++;
	    adjust_diff[ki] = diff_track;
	    	    
    }

    for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
	    if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != exon->strand ) ) continue;
	    if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == exon->strand ) ) continue;
		    
	    uint32_t type = sense_strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;
	    int64_t ue_end_q = exon->end_q - n_backsearch ;
	    for (int64_t ue_end_t = min_ue_end_t; ue_end_t <= max_ue_end_t; ue_end_t++, ue_end_q++) {
		    if (seq->len - ue_end_q -1 < opt->min_anchor) break;
		    bool splice_site_flag = 0; 
		    int64_t intron_start_t = ue_end_t + 1;

		    splice_site_t* known_ss= retrieve_splice_site (opt->splice_site_map, intron_start_t, sense_strand, type);
		    //still test this since the other part might be novel.
		    if (known_ss) {
			    splice_site_flag =1;
		    }
		    else if (jigsaw_is_splice_site_pos (intron_start_t, type, sense_strand, l_pac, pacseq)) {
			    splice_site_flag =1;
		    }

		    if(splice_site_flag) {
			    //init' the anchor seq
			    jigsaw_anchor_seq_t *anchor_seq = (jigsaw_anchor_seq_t*) calloc (1, sizeof(jigsaw_anchor_seq_t) );
			    //get the partner splice site intron pattern
			    ubyte_t *ss = (ubyte_t*) calloc (2, sizeof(ubyte_t) );
			    jigsaw_get_partner_splice_site(type, sense_strand, ss);
							    
			    //get the pos of the other part of the read
			    int64_t de_start_q = ue_end_q + 1;
			    int64_t de_end_q = (int64_t) seq->len -1;
			    anchor_seq->len = anchor_seq->full_len = seq->len - ue_end_q -1 + 2;
			    anchor_seq->seq = (ubyte_t*)calloc(anchor_seq->len, 1);
			    
			    //put the ss into the begining of the anchor
			    int j=0;
			    for (int i = 0; i<2; i++,j++)
			    {
				    anchor_seq->seq[j] = ss[i];
			    }

			    free (ss);
			    // put the other part of the read into anchorseq
			    for (int i = de_start_q; i< de_end_q+1; i++,j++ )
			    {
				    anchor_seq->seq[j] = query_seq [i];
			    }
			    //need to do some init'

			    anchor_seq->rseq = (ubyte_t*)calloc( anchor_seq->len, 1);
			    memcpy(anchor_seq->rseq, anchor_seq->seq, anchor_seq->len);
			    seq_reverse(anchor_seq->len, anchor_seq->rseq, opt->mode & BWA_MODE_COMPREAD);
			    
			    anchor_seq->wid = 0;
			    anchor_seq->offset = 0;
			    anchor_seq->tid = -1; 
			    anchor_seq->qual = 0;
			    //TODO: should copy the qual later
									    
			    // find the hits of this anchor on genome
			    // modified from jigsaw_cal_sa_reg_gap
			    
			    list<jigsaw_word_hit_t*> hits;

			    jigsaw_collect_anchor_hits (bwt, anchor_seq, g_log_n,
					    l_pac, pacseq, ntpac, opt->max_word_occ, opt, &hits);
			    free (anchor_seq->seq);
			    free (anchor_seq->rseq);
			    free (anchor_seq);
			    
			    //check if there is hit
			    if(hits.size() == 0 ) continue;
			    // scan all hits to find optimal junctions
			    // TODO:require perfect match for the anchor for now, can add more things here
			    jigsaw_sort_word_hits_by_strand_pos_t(&hits);

			    int64_t left_bound_t = ue_end_t + opt->min_intron_size;
			    int64_t right_bound_t = ue_end_t + opt->max_intron_size;
		    
			    //fprintf(stderr, "left%u\t%u\n", int(left_bound_t), ue_end_q);	
			    list<jigsaw_word_hit_t*>::iterator iter = hits.begin();
			    //int hit_accepted_flag = 0;
			    for ( ; iter != hits.end(); iter++)
			     {
				jigsaw_word_hit_t *p = *iter;
				int64_t p_pos_t = p->pos_t; 
				if (p->strand == 0 ) {
				    //i think only the seq not the rseq should be mapped
				    if ( p_pos_t <right_bound_t ) {
					if (p_pos_t > left_bound_t ){
					    //create a downstream exon
					    de = (jigsaw_exon_t *)calloc(1, sizeof(jigsaw_exon_t));
					    de->n_mm= de->n_gapo_t = de->n_gape_t = de->n_gapo_q = de->n_gape_q = 0;
					    de->colinear = 1; de->coverage = 1;
					    de->strand = exon->strand;
					    de->start_t = p_pos_t +2;//'cause the splice site pattern
					    
					    //de->end_t = piter->pos_t + anchor_seq->len -1 ;
					    de->end_t = p_pos_t +  seq->len - ue_end_q ;
					    de->start_q = ue_end_q + 1;
					    de->end_q = seq->len -1;
					    de->hits = 0;
					    de->is_first = 0;
					    de->is_last = 1;

					    //create a junction
					    //exon->end_q = ue_end_q;
					    //exon->end_t = ue_end_t;
					    exon->is_last = 0;
					    junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
					    junction->uexon = exon;
					    junction->dexon = de;
					    junction->sense_strand = sense_strand;
					    junction->an_mm = 0;
					    if (ue_end_t - exon->end_t <=0) junction->an_mm = adjust_diff[ exon->end_t-ue_end_t ];
					    junction->n_mm = junction->n_gapo_t = junction->n_gapo_q = 0;
					    junction->n_gape_t = junction->n_gape_q = 0;
					    junction->start_q = ue_end_q;
					    junction->end_q = junction->start_q + 1;
					    junction->start_t = ue_end_t;
					    junction->end_t = de->start_t;
					    junction->logistic_prob = get_splice_score(pacseq, ue_end_t, de->start_t);

					    junctions->push_back(junction);
					    exons->push_back(de);
					    num_junc_found_denovo++;
					}
				    }
				    else break; //exceed the right boundary
				}
				else break; //not on the same strand
			     }

			    //need to free all the hits;

			     for (iter = hits.begin() ; iter != hits.end(); iter++)
			     {
				     jigsaw_word_hit_t *p = *iter;
				     free(p);
			     }
		    }

	    }
    }
    delete[] adjust_diff;
    return num_junc_found_denovo;

}


int jigsaw_locate_junc_one_anchor_with_anno_upstream(jigsaw_exon_t *exon, bwa_seq_t *seq, int64_t min_de_start_t, int64_t max_de_start_t, int n_backsearch, const ubyte_t *pacseq, const gap_opt_t *opt, list<jigsaw_junction_t*> *junctions, list<jigsaw_exon_t*> *exons)
{
    int local_max_diff = (int) opt->max_diff;
    if(local_max_diff>2) local_max_diff = 2;
    //int best_diff = opt->max_diff+1;
    int best_diff = local_max_diff;

    int num_junc_found_in_anno = 0;
    int64_t k;
    jigsaw_junction_t *junction = 0;
    jigsaw_exon_t *ue = 0;
    
     //analog to the two_anchor search
    AlnParam ap = aln_param_bwa;
    //ap.band_width = opt->max_diff;
    ap.band_width = local_max_diff > 1 ? local_max_diff : 1;
    
    //int gap_start_q = 0; //gap starts from beginning of the read
    int gap_len = exon->start_q + n_backsearch + 1;
    
    ubyte_t *query_seq = exon->strand ? seq->rseq : seq->seq;
    
    ubyte_t *gap_seq = query_seq;
    ubyte_t *ref_seq = (ubyte_t *) malloc(gap_len * sizeof(ubyte_t));

    //use an array to store the over counted difference in the boundry of the exon
    uint32_t *adjust_diff =new uint32_t[ n_backsearch +1 ];
    int diff_track = 0;
    for ( int ki =0; ki< n_backsearch + 1; ki++)
    {
	int64_t kt =  exon->start_t + ki;
	int64_t kq = exon->start_q + ki;
	ubyte_t t = get_pacseq_base (pacseq, kt);
	ubyte_t q = query_seq [kq];
	if (q != t) diff_track++;
	adjust_diff[ki] = diff_track;
    }
    //sense_strand: 0 is "+", 1 is "-"
    for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
	    if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != exon->strand ) ) continue;
	    if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == exon->strand ) ) continue;
		    
	    int64_t de_start_q = exon->start_q + n_backsearch;
	    uint32_t type = sense_strand ? SPLICE_DONOR : SPLICE_ACCEPTOR;

	    for (int64_t de_start_t = max_de_start_t; de_start_t >= min_de_start_t ; de_start_t--, de_start_q--) {
		    if (de_start_q  < opt->known_junc_min_anchor) break;
		    int64_t intron_end_t = de_start_t - 1;
		    splice_site_t* known_ss= retrieve_splice_site (opt->splice_site_map, intron_end_t, sense_strand, type);
		    if(!known_ss) continue;
		    //check all the partner sites recorded in splice_site_map
		    int i=gap_len-1;
		    //put the 2nd half of the sequence into ref_seq
		    for (k = max_de_start_t; k >= de_start_t; --k, --i)
			    ref_seq[i] = get_pacseq_base (pacseq,k);
		    for(int partner_i = 0; partner_i < known_ss->n_partners; partner_i++ ){
			    splice_site_t* partner_ss = known_ss->partners[partner_i];
			    int64_t ue_end_t = partner_ss->pos -1;
			    if (ue_end_t - de_start_q  < 0 ) continue;
			    // bug 140119, alignment out of reference boundary
			    int ii=i;
			    //put the other half into ref_seq
			    for (k = ue_end_t; ii>=0; --ii, --k)
				    ref_seq[ii] = get_pacseq_base (pacseq,k);
			    path_t *path = (path_t*)calloc(gap_len * 2, sizeof(path_t));
			    int path_len;
			    aln_global_core(ref_seq, gap_len, gap_seq, gap_len, &ap, path, &path_len);
			    //calculate the number of mismatches and gaps
			    uint32_t n_mm, n_gapo_t, n_gapo_q, n_gape_t, n_gape_q, an_mm;
			    jigsaw_cal_diff (path, path_len, ref_seq, gap_len, gap_seq, gap_len,
			    &n_mm, &n_gapo_t, &n_gapo_q, &n_gape_t, &n_gape_q);
			    free (path);
			    an_mm = 0;
			    if (de_start_t  - exon->start_t >=0) an_mm = adjust_diff[de_start_t - exon->start_t];
			    int diff = n_mm + n_gapo_t + n_gapo_q + n_gape_t + n_gape_q - an_mm;
			    //better junction found
			    if(opt->report_best_only){
				if (diff <= local_max_diff && diff < best_diff) {
					best_diff = diff;
					//create a de first
					if (!ue) ue = (jigsaw_exon_t *)calloc(1, sizeof(jigsaw_exon_t));
					ue->n_mm= ue->n_gapo_t = ue->n_gape_t = ue->n_gapo_q = ue->n_gape_q = 0;
					ue->colinear = 1; ue->coverage = 1;
					ue->strand = exon->strand;
					ue->start_t = ue_end_t - de_start_q +1;
					ue->end_t = ue_end_t; 
					ue->start_q = 0;
					ue->end_q = de_start_q -1; 
					ue->hits = 0;
					ue->is_first = 1;
					ue->is_last = 0;

					//modify the downstream exon
					//exon->start_q = de_start_q;
					//exon->start_t = de_start_t;
					exon->is_first = 0;
					if(!junction) junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
					junction->dexon = exon;
					junction->uexon = ue;
					junction->sense_strand = sense_strand;
					junction->n_mm = n_mm;
					junction->an_mm = an_mm;
					junction->n_gapo_t = n_gapo_t;
					junction->n_gapo_q = n_gapo_q;
					junction->n_gape_t = n_gape_t;
					junction->n_gape_q = n_gape_q;
					junction->start_q = ue->end_q;
					junction->end_q = junction->start_q + 1;
					junction->start_t = ue->end_t;
					junction->end_t = de_start_t;
				    //junction->logistic_prob = get_splice_score(pacseq, ue->end_t, de_start_t);
				    junction->logistic_prob = 1;
					
				}
			    }
			    else{
				//push all the results into the lists
				if(diff <= local_max_diff ){
				    ue = (jigsaw_exon_t *)calloc(1, sizeof(jigsaw_exon_t));
				    ue->n_mm= ue->n_gapo_t = ue->n_gape_t = ue->n_gapo_q = ue->n_gape_q = 0;
				    ue->colinear = 1; ue->coverage = 1;
				    ue->strand = exon->strand;
				    ue->start_t = ue_end_t - de_start_q +1;
				    ue->end_t = ue_end_t; 
				    ue->start_q = 0;
				    ue->end_q = de_start_q -1; 
				    ue->hits = 0;
				    ue->is_first = 1;
				    ue->is_last = 0;

				    //modify the downstream exon
				    //exon->start_q = de_start_q;
				    //exon->start_t = de_start_t;
				    exon->is_first = 0;
				    junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
				    junction->dexon = exon;
				    junction->uexon = ue;
				    junction->sense_strand = sense_strand;
				    junction->n_mm = n_mm;
				    junction->an_mm = an_mm;
				    junction->n_gapo_t = n_gapo_t;
				    junction->n_gapo_q = n_gapo_q;
				    junction->n_gape_t = n_gape_t;
				    junction->n_gape_q = n_gape_q;
				    junction->start_q = ue->end_q;
				    junction->end_q = junction->start_q + 1;
				    junction->start_t = ue->end_t;
				    junction->end_t = de_start_t;
				    //junction->logistic_prob = get_splice_score(pacseq, ue->end_t, de_start_t);
				    junction->logistic_prob = 1;

				    junctions->push_back(junction);
				    exons->push_back(ue);
				    num_junc_found_in_anno++;
				}
				    

			    }
				
			     
		    }
			    
	    }
    }
    free(ref_seq);
    if(opt->report_best_only && junction){
      junctions->push_back(junction);
      exons->push_back(ue);
      num_junc_found_in_anno++;
    }
    delete adjust_diff;
    return(num_junc_found_in_anno);
		

}

int jigsaw_locate_junc_one_anchor_denovo_upstream(jigsaw_exon_t *exon, bwa_seq_t *seq, int64_t min_de_start_t, int64_t max_de_start_t, int n_backsearch, int64_t l_pac,const ubyte_t *pacseq, const ubyte_t *ntpac, bwt_t *const bwt[2], const int *g_log_n, const gap_opt_t *opt, list<jigsaw_junction_t*> *junctions, list<jigsaw_exon_t*> *exons)
{
	if (! (exon->start_q + opt->max_overhang >= opt->min_anchor && exon->start_q < 2 * seq->word_size)) return 0;
    int num_junc_found_denovo = 0;
    ubyte_t *query_seq = exon->strand ? seq->rseq : seq->seq;
    jigsaw_exon_t *ue = 0;
    jigsaw_junction_t *junction = 0;

    uint32_t *adjust_diff =new uint32_t[ n_backsearch +1 ];
    int diff_track = 0;
    for ( int ki =0; ki< n_backsearch + 1; ki++)
    {
	    int64_t kt =  exon->start_t + ki;
	    int64_t kq = exon->start_q + ki;
	    ubyte_t t = get_pacseq_base (pacseq, kt);
	    ubyte_t q = query_seq [kq];
	    if (q != t) diff_track++;
	    adjust_diff[ki] = diff_track;
    }


    for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
	if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != exon->strand ) ) continue;
	if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == exon->strand ) ) continue;
		    
        uint32_t type = sense_strand ? SPLICE_DONOR : SPLICE_ACCEPTOR; 
	int64_t de_start_q = exon->start_q + n_backsearch ;
	for (int64_t de_start_t = max_de_start_t; de_start_t >= min_de_start_t ; de_start_t--, de_start_q--) {
	    if (de_start_q < opt->min_anchor) break;
	    bool splice_site_flag = 0;
	    int64_t intron_end_t = de_start_t - 1;

	    splice_site_t* known_ss= retrieve_splice_site (opt->splice_site_map, intron_end_t, sense_strand, type);
	    if (known_ss) {
		    splice_site_flag =1;
	    }
	    else if (jigsaw_is_splice_site_pos (intron_end_t, type, sense_strand, l_pac, pacseq)) {
		    splice_site_flag =1;
	    }
	    if(splice_site_flag) {
		    jigsaw_anchor_seq_t *anchor_seq = (jigsaw_anchor_seq_t*) calloc (1, sizeof(jigsaw_anchor_seq_t) );                  
		    //get the partner splice site intron pattern
		    //ubyte_t *ss = jigsaw_get_partner_splice_site (type, sense_strand);
		    ubyte_t *ss = (ubyte_t*) calloc (2, sizeof(ubyte_t) );
		    jigsaw_get_partner_splice_site(type, sense_strand, ss);
		    
		    //get the pos of the other part of the read
		    int64_t ue_end_q = de_start_q -1;
		    int64_t ue_start_q = 0; //start of the read
		    anchor_seq->len = anchor_seq->full_len = de_start_q  +2;
		    anchor_seq->seq = (ubyte_t*)calloc(anchor_seq->len, 1);

		    int j=0;
		    // put the other part of the read into anchorseq
		    for (int i = ue_start_q; i< ue_end_q+1; i++,j++ )
		    {
			    anchor_seq->seq[j] = query_seq [i];
		    }
		    
		    //put the ss into the end of the anchor
		    for (int i = 0; i<2; i++,j++)
		    {
			    anchor_seq->seq[j] = ss[i];
		    }
		    free (ss);
		    anchor_seq->rseq = (ubyte_t*)calloc( anchor_seq->len, 1);
		    memcpy(anchor_seq->rseq, anchor_seq->seq, anchor_seq->len);
		    seq_reverse(anchor_seq->len, anchor_seq->rseq, opt->mode & BWA_MODE_COMPREAD);

		    anchor_seq->wid = 0;
		    anchor_seq->offset = 0;
		    anchor_seq->tid = -1;
		    anchor_seq->qual = 0;
		    //TODO: should copy the qual later
		    list<jigsaw_word_hit_t*> hits;
		    jigsaw_collect_anchor_hits (bwt, anchor_seq, g_log_n,
						    l_pac, pacseq, ntpac, opt->max_word_occ, opt, &hits);
		    //check if there is hit
		    free (anchor_seq->seq);
		    free (anchor_seq->rseq);
		    free (anchor_seq);
		    if(hits.size() == 0 ) continue;
		    // scan all hits to find optimal junctions
		    // sort, then only look at the nearest match
		    // TODO:require perfect match for the anchor for now, can add more things here
		    jigsaw_sort_word_hits_by_strand_pos_t(&hits);
		    int64_t right_bound_t = de_start_t - opt->min_intron_size - de_start_q;
		    int64_t left_bound_t = de_start_t - opt->max_intron_size - de_start_q;
							    
		    list<jigsaw_word_hit_t*>::iterator iter = hits.end();
		    for ( ; iter != hits.begin(); )
		    {
			iter--;
			//go backwards this time
			jigsaw_word_hit_t *p = *iter;
			int64_t p_pos_t = p->pos_t;
			if ( p->strand == 0){
			     if(p_pos_t > left_bound_t){
				 if( p_pos_t <right_bound_t ) {
				     ue = (jigsaw_exon_t *)calloc(1, sizeof(jigsaw_exon_t));
				     ue->n_mm= ue->n_gapo_t = ue->n_gape_t = ue->n_gapo_q = ue->n_gape_q = 0;
				     ue->colinear = 1; ue->coverage = 1;
				     ue->strand = exon->strand;
				     ue->start_t = p_pos_t ;
				     //ue->end_t = piter->pos_t + anchor_seq->len -1 -2 ;
				     ue->end_t = p_pos_t + de_start_q -1 ;

				     ue->start_q = 0;
				     ue->end_q = ue_end_q;
				     ue->hits = 0;
				     ue->is_first = 1;
				     ue->is_last = 0;
				     
				    exon->is_first = 0;
				    junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
				    junction->dexon = exon;
				    junction->uexon = ue;
				    junction->sense_strand = sense_strand;
				    junction->n_mm = junction->n_gapo_t = junction->n_gapo_q = 0;
				    junction->an_mm = 0;
				    if (de_start_t  - exon->start_t >=0) junction->an_mm = adjust_diff[de_start_t - exon->start_t];
				    junction->n_gape_t = junction->n_gape_q = 0;
				    junction->start_q = ue->end_q;
				    junction->end_q = junction->start_q + 1;
				    junction->start_t = ue->end_t;
				    junction->end_t = de_start_t;
				    junction->logistic_prob = get_splice_score(pacseq, ue->end_t, de_start_t);
				    
				    junctions->push_back(junction);
				    exons->push_back(ue);
				    num_junc_found_denovo++;
				 }
			     }
			     else break;//exceed the right boundary
			}
		    }
		    for (iter = hits.begin() ; iter != hits.end(); iter++)
		     {
			     jigsaw_word_hit_t *p = *iter;
			     free(p);
		     }

	    }
		    
	}
    }
    delete[] adjust_diff;

    return num_junc_found_denovo;	

}
/*
   find a splice site 
   -> retrieve the rest part of the read
   -> plus the partner splicepattern 
   -> align with BWT -> output results
*/

void jigsaw_locate_junc_one_anchor (jigsaw_exon_t* exon, bwa_seq_t *seq, int direction, const gap_opt_t *opt,	int64_t l_pac,const ubyte_t *pacseq, const ubyte_t *ntpac, bwt_t *const bwt[2], const int *g_log_n, list<jigsaw_junction_t*> *junctions, list<jigsaw_exon_t*> *exons)
{
	int n_backsearch = opt->max_overhang;
	
	if (direction == 0){
		//direction = 0, downstream intron
		int64_t min_ue_end_t = exon->end_t - n_backsearch;
		int64_t max_ue_end_t = exon->end_t;
		
		//look for 5' splice site downstream
		int num_junc_found_in_anno = 0, num_junc_found_denovo = 0;
		//if a junction database is provided, search with it first
		if(opt->splice_site_map)
		{
		    num_junc_found_in_anno = jigsaw_locate_junc_one_anchor_with_anno_downstream(exon, seq, min_ue_end_t, max_ue_end_t, n_backsearch, l_pac, pacseq, opt, junctions, exons);
		}
		//continue the following denovo search

		if (! opt->non_denovo_search) num_junc_found_denovo = jigsaw_locate_junc_one_anchor_denovo_downstream(exon, seq, min_ue_end_t, max_ue_end_t,  n_backsearch,  l_pac, pacseq, ntpac,  bwt, g_log_n, opt, junctions, exons);
	}
	else {

		//the other direction, everything reverse
		//look for upstream exon
		int64_t max_de_start_t = exon->start_t + n_backsearch;
		int64_t min_de_start_t = exon->start_t;
		
		int num_junc_found_in_anno = 0, num_junc_found_denovo = 0;
		//if a junction database is provided, search with it first
		if(opt->splice_site_map){
		    num_junc_found_in_anno = jigsaw_locate_junc_one_anchor_with_anno_upstream(exon, seq, min_de_start_t, max_de_start_t, n_backsearch, pacseq, opt, junctions, exons);
		}
		if (! opt->non_denovo_search)	num_junc_found_denovo = jigsaw_locate_junc_one_anchor_denovo_upstream(exon, seq, min_de_start_t, max_de_start_t,  n_backsearch,  l_pac, pacseq, ntpac,  bwt, g_log_n, opt, junctions, exons);
		
	}
	return;
}

/*
 * this is the old single anchor search function. 
 * find exon junctions using one anchor (limited to known junctions)
 * direction: look for upstream (1) or downstream (0) exon
 * parter: the partner exon found
 * return the new junction
 */
/*
jigsaw_junction_t* _jigsaw_locate_junc_one_anchor (jigsaw_exon_t* exon, bwa_seq_t *seq, int direction,
		splice_site_map_t* splice_site_map, int max_intron_size, int min_intron_size, int max_overhang, int min_anchor,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac, jigsaw_exon_t* partner)
{
	int n_backsearch = max_overhang;

	//int minOverlap = 4;
	bool stop_if_neg_score = false;
	jigsaw_junction_t* junction = 0;

	if (direction == 0) {
		//direction = 0, downstream intron

		int64_t min_ue_end_t = exon->end_t - n_backsearch;

		int64_t max_ue_end_t = exon->end_t;
		//int32_t strand = exon->strand;
		int64_t ue_end_t;

		jigsaw_exon_t *de = 0;
		//int max_de_len = 0;
		int min_de_mm = seq->len; //TODO: to be replaced

		//sense_strand: 0 is "+", 1 is "-"
		for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
			uint32_t type = sense_strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;

			//make sure there is at least one splice site exists. otherwise we don't need to enumerate partner splice sites
			int c = 0;
			for (ue_end_t = min_ue_end_t; ue_end_t <= max_ue_end_t; ++ue_end_t) {
				int64_t intron_start_t = ue_end_t + 1;

				splice_site_t* known_ss= retrieve_splice_site (splice_site_map, intron_start_t, sense_strand, type);
				if (known_ss) {
					c++; break;
				}

				if (jigsaw_is_splice_site_pos (intron_start_t, type, sense_strand, l_pac, pacseq)) {
					c++; break;
				}
			}

			if (c < 1) continue;


			//enumerate partner splice sites de novo
			list <splice_site_t *> partner_splice_sites;

			uint32_t partner_type = sense_strand ? SPLICE_DONOR : SPLICE_ACCEPTOR;

			int64_t enum_start_t = min_ue_end_t > 0 ? min_ue_end_t+1 : 1;
			int64_t enum_end_t = (max_ue_end_t+max_intron_size) < l_pac ? max_ue_end_t+max_intron_size: l_pac;

			jigsaw_enum_splice_site_denovo (enum_start_t, enum_end_t, sense_strand, partner_type, direction,
					l_pac, pacseq, &partner_splice_sites);


			//now ready to do the search
			for (ue_end_t = min_ue_end_t; ue_end_t <= max_ue_end_t; ++ue_end_t) {

				int64_t intron_start_t = ue_end_t + 1;

				//list <splice_site_t *> local_partner_splice_sites;

				//record the original end
				list <splice_site_t *>::iterator old_end = partner_splice_sites.end ();

				//add novel partner splice sites only if this is a valid splice site
				//if (jigsaw_is_splice_site_pos (intron_start_t, type, sense_strand, l_pac, pacseq))
				//	local_partner_splice_sites.assign (partner_splice_sites.begin(), partner_splice_sites.end());

				int c = 0;
				//add known partner splice sites
				if (splice_site_map)
				{
					splice_site_t* ss= retrieve_splice_site (splice_site_map, intron_start_t, sense_strand, type);
					if (ss)
					{
						c++;
						splice_site_t **partners = ss->partners;
						for (int i = 0; i < ss->n_partners; ++i) {
							partner_splice_sites.push_back (partners[i]);
						}
					}
				}

				if (jigsaw_is_splice_site_pos (intron_start_t, type, sense_strand, l_pac, pacseq)) c++;

				if (c < 1) continue; // make sure this is a legal splice site

				for (list <splice_site_t *>::iterator iter = partner_splice_sites.begin(); iter!= partner_splice_sites.end(); ++iter)
				{
					splice_site_t *p = *iter;
					int64_t intron_end_t = p->pos;

					//make sure the intron is downstream of the exon and has acceptable size
					if (!(intron_end_t - intron_start_t + 1 >= min_intron_size && intron_end_t - intron_start_t + 1 <= max_intron_size)) continue;

					//create a new downstream exon
					if (!de) de = (jigsaw_exon_t *)calloc(1, sizeof(jigsaw_exon_t));
					de->n_mm= de->n_gapo_t = de->n_gape_t = de->n_gapo_q = de->n_gape_q = 0;
					de->colinear = 1; de->coverage = 1;
					de->strand = exon->strand;
					de->start_t = de->end_t = intron_end_t + 1;
					de->start_q = de->end_q = exon->end_q - n_backsearch + ue_end_t - min_ue_end_t + 1;
					de->hits = 0; //no hits
					de->is_first = 0; de->is_last = 1;

					//extend the exon on the right
					jigsaw_extend_exon_inexact (de, seq, 0, l_pac, pacseq, ntpac, stop_if_neg_score); //extend the exon on the right

					int de_len = de->end_q - de->start_q + 1; // - de->n_mm; //corrected length

					//create a junction
					//TODO: right now junctions are ranked according to the number of nucleotides on the other end
					//need to consider the number of mismatches as well
					int de_mm = de->n_mm + seq->len - 1- de->end_q;
					//bool better_de = de->n_mm < min_de_mm || (de->n_mm == min_de_mm && de_len > max_de_len);

					if (de_len >= min_anchor && de_mm < min_de_mm)
					{
						min_de_mm = de_mm;
						//max_de_len = de_len;

						*partner = *de;
						exon->is_last = 0;

						if (junction) free (junction);
						junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
						junction->uexon = exon;
						junction->dexon = partner;
						junction->strand = sense_strand;
						junction->n_mm = junction->n_gapo_t = junction->n_gapo_q = 0;
						junction->n_gape_t = junction->n_gape_q = 0;
						junction->start_q = exon->end_q - n_backsearch + ue_end_t - min_ue_end_t;
						junction->end_q = junction->start_q + 1;
						junction->start_t = ue_end_t; junction->end_t = de->start_t;
					}

					//free (p);

					if (de->end_q == seq->len - 1 && de->n_mm == 0)
					{
						//break;
						partner_splice_sites.erase (old_end, partner_splice_sites.end ());
						for (list <splice_site_t *>::iterator iter = partner_splice_sites.begin();
							iter!= partner_splice_sites.end(); ++iter)
						{
							splice_site_t *p = *iter;
							free (p);
						}
						if (de) free (de);
						return junction;
					}
					//TODO: if the best possible partner is already found, stop now
					//need to pair the new exon with other downstream exons
				} //loop of partners


				partner_splice_sites.erase (old_end, partner_splice_sites.end ());

			}//loop of ue_end

			for (list <splice_site_t *>::iterator iter = partner_splice_sites.begin();
						iter!= partner_splice_sites.end(); ++iter)
			{
				splice_site_t *p = *iter;
				free (p);
			}

		}//loop of strand
		if (de) free (de);
	}
	else {
		//direction = 1, upstream intron

		int64_t min_de_start_t = exon->start_t;
		int64_t max_de_start_t = exon->start_t + n_backsearch;
		//int32_t strand = exon->strand;
		int64_t de_start_t;

		jigsaw_exon_t *ue = 0;
		//int max_ue_len = 0;
		int min_ue_mm = seq->len; //TODO: to be replaced

		for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
			uint32_t type = sense_strand ? SPLICE_DONOR : SPLICE_ACCEPTOR;

			//make sure there is at least one splice site exists. otherwise we don't need to enumerate partner splice sites
			int c = 0;
			for (de_start_t = min_de_start_t; de_start_t <= max_de_start_t; ++de_start_t) {
				int64_t intron_end_t = de_start_t - 1;

				splice_site_t* known_ss= retrieve_splice_site (splice_site_map, intron_end_t, sense_strand, type);
				if (known_ss) {
					c++; break;
				}

				if (jigsaw_is_splice_site_pos (intron_end_t, type, sense_strand, l_pac, pacseq)) {
					c++;break;
				}
			}

			if (c < 1) continue;

			//enumerate partner splice sites de novo
			list <splice_site_t *> partner_splice_sites;

			uint32_t partner_type = sense_strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;
			int64_t enum_start_t = max_de_start_t > max_intron_size ? max_de_start_t - max_intron_size: 1;
			int64_t enum_end_t = (max_de_start_t - 1) < l_pac ? max_de_start_t - 1 : l_pac;
			
			jigsaw_enum_splice_site_denovo ( enum_start_t, enum_end_t, sense_strand, partner_type, direction,
					l_pac, pacseq, &partner_splice_sites);

			//now ready to search
			for (de_start_t = min_de_start_t; de_start_t <= max_de_start_t; ++de_start_t) {
				int64_t intron_end_t = de_start_t - 1;

				//record the original end
				list <splice_site_t *>::iterator old_end = partner_splice_sites.end ();


				//list <splice_site_t *> local_partner_splice_sites;

				//add novel partner splice sites only if this is a valid splice site
				//if (jigsaw_is_splice_site_pos (intron_end_t, type, sense_strand, l_pac, pacseq))
				//	local_partner_splice_sites.assign (partner_splice_sites.begin(), partner_splice_sites.end());

				int c = 0;

				//add known splice sites
				if (splice_site_map)
				{
					splice_site_t* ss = retrieve_splice_site (splice_site_map, intron_end_t, sense_strand, type);
					if (ss)
					{
						c++;
						splice_site_t **partners = ss->partners;
						for (int i = 0; i < ss->n_partners; ++i)
						{
							partner_splice_sites.push_back (partners[i]);
						}
					}
				}

				if (jigsaw_is_splice_site_pos (intron_end_t, type, sense_strand, l_pac, pacseq)) c++;

				if (c < 1) continue;


				for (list <splice_site_t *>::iterator iter = partner_splice_sites.begin(); iter!= partner_splice_sites.end(); ++iter)
				{
					splice_site_t *p = *iter;
					int64_t intron_start_t = p->pos;

					//fprintf (stderr, "intron start=%llu, intron end=%llu\n", intron_start_t, intron_end_t);


					//make sure the intron is downstream of the exon and has acceptable size
					if (!(intron_end_t - intron_start_t + 1 >= min_intron_size && intron_end_t - intron_start_t + 1 <= max_intron_size)) continue;

					//create a new upstream exon
					if (!ue) ue = (jigsaw_exon_t *)calloc(1, sizeof(jigsaw_exon_t));
					ue->n_mm= ue->n_gapo_t = ue->n_gape_t = ue->n_gapo_q = ue->n_gape_q = 0;
					ue->colinear = 1; ue->coverage = 1;
					ue->strand = exon->strand;
					ue->start_t = ue->end_t = intron_start_t - 1;
					ue->start_q = ue->end_q = exon->start_q + de_start_t - min_de_start_t - 1;
					ue->hits = 0; //no hits
					ue->is_first = 1; ue->is_last = 0;
					//extend the exon on the left
					jigsaw_extend_exon_inexact (ue, seq, 1, l_pac, pacseq, ntpac, stop_if_neg_score);
					int ue_len = ue->end_q - ue->start_q + 1; // - ue->n_mm;

					int ue_mm = ue->n_mm + ue->start_q;
					//bool better_ue = ue->n_mm < min_ue_mm || (de->n_mm == min_de_mm && de_len > max_de_len);


					//create a junction
					if (ue_len >= min_anchor && ue_mm < min_ue_mm) {
						//max_ue_len = ue_len;
						min_ue_mm = ue_mm;

						*partner = *ue;

						if (junction) free (junction);
						junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
						junction->uexon = partner;
						junction->dexon = exon;
						junction->strand = sense_strand;
						junction->n_mm = junction->n_gapo_t = junction->n_gapo_q = 0;
						junction->n_gape_t = junction->n_gape_q = 0;
						junction->start_q = ue->end_q;
						junction->end_q = junction->start_q + 1;
						junction->start_t = ue->end_t; junction->end_t = de_start_t;
					}

					//free (p);

					if (ue->start_q == 0 && ue->n_mm == 0)
					{
						//TODO: if the best possible partner is already found, stop now
						//need to delete the known splice site from the local list and release the memory

						partner_splice_sites.erase (old_end, partner_splice_sites.end ());
						for (list <splice_site_t *>::iterator iter = partner_splice_sites.begin();
							iter!= partner_splice_sites.end(); ++iter)
						{
							splice_site_t *p = *iter;
							free (p);
						}
						if (ue) free (ue);
						return junction;
					}
					//need to pair the new exon with other upstream exons
				} //loop of partners

				partner_splice_sites.erase (old_end, partner_splice_sites.end ());
			}//loop of de_start

			for (list <splice_site_t *>::iterator iter = partner_splice_sites.begin();
				iter!= partner_splice_sites.end(); ++iter)
			{
				splice_site_t *p = *iter;
				free (p);
			}
		}//loop of strand
		if (ue) free (ue);
	}
	return junction;
}
*/

//when a splice site annotation is provided,
//return the number of junctions found
int jigsaw_locate_junc_two_anchors_with_anno (jigsaw_exon_t *ue, jigsaw_exon_t *de, int64_t min_ue_end_t, int64_t max_ue_end_t, int n_backsearch, uint32_t seq_len, int gap_len, ubyte_t *seq, ubyte_t *ref_seq,  const ubyte_t *pacseq, const AlnParam *ap, const gap_opt_t *opt, list<jigsaw_junction_t*> *junctions)

{
	ubyte_t *gap_seq = seq + ue->end_q - n_backsearch + 1;
	ubyte_t *query_seq = seq;
	int local_max_diff = (int) opt->max_diff;
	if(local_max_diff >2 ) local_max_diff = 2;
        //int best_diff = opt->max_diff+1;
	int best_diff = local_max_diff + 1;
	int i;
	int64_t k;
	jigsaw_junction_t* p = 0;
	int num_junc_found_in_anno = 0;

	uint32_t *ue_adjust_diff =new uint32_t[ n_backsearch +1 ];
	int diff_track = 0;
	for ( int ki =0; ki< n_backsearch + 1; ki++) {
		int64_t kt = ue->end_t -ki;
		int64_t kq = ue->end_q - ki;
		ubyte_t t = get_pacseq_base (pacseq, kt);
		ubyte_t q = query_seq [kq];
		if (q != t) diff_track++;
		ue_adjust_diff[ki] = diff_track;
	}
	uint32_t *de_adjust_diff = new uint32_t [  n_backsearch +1 ];
	diff_track = 0;
	for ( int ki =0; ki< n_backsearch + 1; ki++) {
		int64_t kt = de->start_t + ki;
		int64_t kq = de->start_q + ki;
		ubyte_t t = get_pacseq_base (pacseq, kt);
		ubyte_t q = query_seq [kq];
		if (q != t) diff_track++;
		de_adjust_diff[ki] = diff_track;
	}

	for (int64_t ue_end_t = min_ue_end_t; ue_end_t <= max_ue_end_t; ++ue_end_t) {
		//remove positions which are too close to the boundary 
                int breakpoint_pos_q_tmp = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
                if(breakpoint_pos_q_tmp +1 < opt->known_junc_min_anchor || (int (seq_len) - breakpoint_pos_q_tmp) < opt->known_junc_min_anchor ) continue;
	       int64_t de_start_t = (de->start_t + n_backsearch) - (gap_len - (ue_end_t - min_ue_end_t));
	       for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
		   if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != ue->strand )  ) continue;
		   if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == ue->strand ) ) continue;
		   
		  uint32_t type = sense_strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;
		  splice_site_t* known_ss = retrieve_splice_site (opt->splice_site_map, ue_end_t +1, sense_strand, type);
		  if (!known_ss) continue;
		  bool match_an_intron = 0; 
		  for(int partner_i = 0; partner_i < known_ss->n_partners; partner_i++ )
		  {
			  splice_site_t* partner_ss = known_ss->partners[partner_i];
			  if(partner_ss->pos == (unsigned)de_start_t-1) {
				  match_an_intron = 1;
				  break;
			  }
		  }

		  if(match_an_intron){
			  //perform banded alignment near candidate splice sites
			  i = 0;
			  for (k = min_ue_end_t + 1; k <= ue_end_t; ++k, ++i)
				  ref_seq[i] = get_pacseq_base (pacseq,k);
			  for (k = de_start_t; k < de->start_t + n_backsearch; ++k, ++i)
				  ref_seq[i] = get_pacseq_base (pacseq,k);
			  path_t *path = (path_t*)calloc(gap_len * 2, sizeof(path_t));
			  int path_len;
			  aln_global_core(ref_seq, gap_len, gap_seq, gap_len, ap, path, &path_len);
			  //calculate the number of mismatches and gaps
			  uint32_t n_mm, n_gapo_t, n_gapo_q, n_gape_t, n_gape_q;
			  jigsaw_cal_diff (path, path_len, ref_seq, gap_len, gap_seq, gap_len,
			   &n_mm, &n_gapo_t, &n_gapo_q, &n_gape_t, &n_gape_q);
			  free (path);
			  uint32_t an_mm = 0;
			  //if (ue_end_t - ue->end_t <=0)  an_mm +=  ue_adjust_diff[ ue->end_t-ue_end_t ];
			  //if (de_start_t  - de->start_t >=0)   an_mm += de_adjust_diff[de_start_t - de->start_t];
			  an_mm += ue_adjust_diff[ n_backsearch -1 ];
			  an_mm += de_adjust_diff[ n_backsearch -1 ];

			  int diff = n_mm + n_gapo_t + n_gapo_q + n_gape_t + n_gape_q - an_mm;

			  if(opt->report_best_only){  
			      //better junction found
			      if (diff <= local_max_diff && diff < best_diff) {
					if (!p) p = (jigsaw_junction_t *) calloc (1, sizeof(jigsaw_junction_t));
					p->sense_strand = sense_strand; //strand of the splice site
					p->uexon = ue; p->dexon = de;
					p->start_t = ue_end_t; p->end_t = de_start_t;
					p->start_q = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
					p->end_q = p->start_q + 1;

					p->n_mm = n_mm;
					p->an_mm = an_mm;
					p->n_gapo_t = n_gapo_t; p->n_gape_t = n_gape_t;
					p->n_gapo_q = n_gapo_q; p->n_gape_q = n_gape_q;
					
					//p->logistic_prob = get_splice_score(pacseq, ue_end_t, de_start_t);
					p->logistic_prob = 1;
					ue->is_last = 0; de->is_first = 0;
					
					best_diff = diff;

			      }
			 }
			 else {
			    //push all results into the junctions
			    if (diff <= local_max_diff ){
				p = (jigsaw_junction_t *) calloc (1, sizeof(jigsaw_junction_t));
				
				p->sense_strand = sense_strand;
				p->uexon = ue; p->dexon = de;
				//a number of junctions will share these exons....for now
				p->start_t = ue_end_t; p->end_t = de_start_t;
				p->start_q = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
				p->end_q = p->start_q + 1;
				p->n_mm = n_mm;
				p->an_mm = an_mm;
				p->n_gapo_t = n_gapo_t; p->n_gape_t = n_gape_t;
				p->n_gapo_q = n_gapo_q; p->n_gape_q = n_gape_q;
				
				ue->is_last = 0; de->is_first = 0;	
				//p->logistic_prob = get_splice_score(pacseq, ue_end_t, de_start_t);
				p->logistic_prob = 1;
							
				junctions->push_back(p);  
				num_junc_found_in_anno++;
			    }
			 }
		  }
		  else continue;
	  
	       }

	    }
	    
	    if(opt->report_best_only && p){
		junctions->push_back(p);
		num_junc_found_in_anno++;
		
	    }
	delete[] ue_adjust_diff;
	delete[] de_adjust_diff;
	return num_junc_found_in_anno;    
}

//this function will do denovo search of SINGLE exon between two anchors
//will return the number of junctions found denovo
int jigsaw_locate_junc_two_anchors_denovo (jigsaw_exon_t *ue, jigsaw_exon_t *de, int64_t min_ue_end_t, int64_t max_ue_end_t, int n_backsearch, uint32_t seq_len, int gap_len, ubyte_t *seq, ubyte_t *ref_seq,  const ubyte_t *pacseq,const AlnParam *ap, const gap_opt_t *opt, list<jigsaw_junction_t*> *junctions)
{
	ubyte_t *gap_seq = seq + ue->end_q - n_backsearch + 1;
	ubyte_t *query_seq = seq;
        //int local_max_diff = int (opt->max_diff/2 +0.5);
	int local_max_diff = int (opt->max_diff);
	if( local_max_diff > 2 ) local_max_diff = 2;
        //int best_diff = opt->max_diff+1;
	int best_diff = local_max_diff +1;
	uint32_t sense_strand; //strand of splice sites
	int num_junc_found_denovo = 0;
	int64_t istart_t, iend_t,k ;
	int i;
	jigsaw_junction_t* p = 0;
        uint32_t *ue_adjust_diff =new uint32_t[ n_backsearch +1 ];
        int diff_track = 0;
        for ( int ki =0; ki< n_backsearch + 1; ki++) {
                int64_t kt = ue->end_t -ki;
                int64_t kq = ue->end_q - ki;
                ubyte_t t = get_pacseq_base (pacseq, kt);
                ubyte_t q = query_seq [kq];
                if (q != t) diff_track++;
                ue_adjust_diff[ki] = diff_track;
        }
        uint32_t *de_adjust_diff = new uint32_t [  n_backsearch +1 ];
        diff_track = 0;
        for ( int ki =0; ki< n_backsearch + 1; ki++) {
                int64_t kt = de->start_t + ki;
                int64_t kq = de->start_q + ki;
                ubyte_t t = get_pacseq_base (pacseq, kt);
                ubyte_t q = query_seq [kq];
                if (q != t) diff_track++;
                de_adjust_diff[ki] = diff_track;
        }

	
        for (int64_t ue_end_t = min_ue_end_t; ue_end_t <= max_ue_end_t; ++ue_end_t) {
		int64_t de_start_t = (de->start_t + n_backsearch) - (gap_len - (ue_end_t - min_ue_end_t));
		if (de_start_t - ue_end_t < opt->min_intron_size) continue;
		istart_t = ue_end_t + 1; iend_t = de_start_t - 1;
		//sometimes the upstream part or downstream part is too small, for example 85,5 for 90nt read, try to deliminate them here
		int breakpoint_pos_q_tmp = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
		if(breakpoint_pos_q_tmp +1 < opt->min_anchor || (int (seq_len) - breakpoint_pos_q_tmp) < opt->min_anchor ) continue;

		//check splice sites
		//TODO: only GT/AG splice sites now
		if (get_pacseq_base (pacseq,istart_t) == BASE_G && get_pacseq_base (pacseq, istart_t +1) == BASE_T
			&& get_pacseq_base (pacseq, iend_t-1) == BASE_A && get_pacseq_base (pacseq, iend_t) == BASE_G  )
			sense_strand = 0;
		else if (get_pacseq_base (pacseq,istart_t) == BASE_C && get_pacseq_base (pacseq, istart_t +1) == BASE_T
			&& get_pacseq_base (pacseq, iend_t-1) == BASE_A && get_pacseq_base (pacseq, iend_t) == BASE_C  )
			sense_strand = 1;
		else continue;
		
		if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != ue->strand )  ) continue;
                if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == ue->strand ) ) continue;
		   
		//perform banded alignment near candidate splice sites
		i = 0;
		for (k = min_ue_end_t + 1; k <= ue_end_t; ++k, ++i)
			ref_seq[i] = get_pacseq_base (pacseq,k);

		for (k = de_start_t; k < de->start_t + n_backsearch; ++k, ++i)
			ref_seq[i] = get_pacseq_base (pacseq,k);

		path_t *path = (path_t*)calloc(gap_len * 2, sizeof(path_t));
		int path_len;
		aln_global_core(ref_seq, gap_len, gap_seq, gap_len, ap, path, &path_len);

		//calculate the number of mismatches and gaps
		uint32_t n_mm, n_gapo_t, n_gapo_q, n_gape_t, n_gape_q;
		jigsaw_cal_diff (path, path_len, ref_seq, gap_len, gap_seq, gap_len,
			&n_mm, &n_gapo_t, &n_gapo_q, &n_gape_t, &n_gape_q);
		free (path);
		uint32_t an_mm = 0;
		//if (ue_end_t - ue->end_t <=0)  an_mm +=  ue_adjust_diff[ ue->end_t-ue_end_t ];
		//if (de_start_t  - de->start_t >=0)   an_mm += de_adjust_diff[de_start_t - de->start_t];

		//due to the local band alignment here, the an_mm is computed differently.
		an_mm += ue_adjust_diff[ n_backsearch -1 ]; // note the -1
		 an_mm += de_adjust_diff[ n_backsearch -1];

		int diff = n_mm + n_gapo_t + n_gapo_q + n_gape_t + n_gape_q - an_mm;

		//better junction found
		if(opt->report_best_only){
		    if (diff <= local_max_diff && diff < best_diff) {
			    if (!p) p = (jigsaw_junction_t *) calloc (1, sizeof(jigsaw_junction_t));
			    //only allocate once in this case
			    p->sense_strand = sense_strand;
			    p->uexon = ue; p->dexon = de;
			    p->start_t = ue_end_t; p->end_t = de_start_t;
			    p->start_q = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
			    p->end_q = p->start_q + 1;

			    p->n_mm = n_mm;
			    p->an_mm = an_mm;
			    p->n_gapo_t = n_gapo_t; p->n_gape_t = n_gape_t;
			    p->n_gapo_q = n_gapo_q; p->n_gape_q = n_gape_q;
			    p->logistic_prob = get_splice_score(pacseq, ue_end_t, de_start_t);

			    best_diff = diff;
			    ue->is_last = 0; de->is_first = 0;
		    }
		}
		else {
		    if (diff <= local_max_diff){
			p = (jigsaw_junction_t *) calloc (1, sizeof(jigsaw_junction_t));
			p->sense_strand = sense_strand;
			p->uexon = ue; p->dexon = de;
			p->start_t = ue_end_t; p->end_t = de_start_t;
			p->start_q = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
			p->end_q = p->start_q + 1;
			p->n_mm = n_mm;
			p->an_mm = an_mm;
			p->n_gapo_t = n_gapo_t; p->n_gape_t = n_gape_t;
			p->n_gapo_q = n_gapo_q; p->n_gape_q = n_gape_q;
			ue->is_last = 0; de->is_first = 0;
			p->logistic_prob = get_splice_score(pacseq, ue_end_t, de_start_t);
			
			junctions->push_back(p);
			
			num_junc_found_denovo++;
		    }
		}
	}
	//this is the end of denovo search
	if(opt->report_best_only && p){
	    junctions->push_back(p);
	    num_junc_found_denovo++;
	}
	delete[] ue_adjust_diff;
	delete[] de_adjust_diff;
	return num_junc_found_denovo;
} 


int jigsaw_locate_junc_two_anchors_inner_exon_with_anno (jigsaw_exon_t *ue, jigsaw_exon_t *de, int64_t min_ue_end_t, int64_t max_ue_end_t, int n_backsearch,uint32_t seq_len,  ubyte_t *seq, const ubyte_t *pacseq,const ubyte_t *ntpac, int64_t l_pac, bwt_t *const bwt[2], const int *g_log_n, const gap_opt_t *opt, list<jigsaw_junction_t*> *junctions, list<jigsaw_exon_t*> *exons)
{
    ubyte_t *query_seq = seq; 
    int64_t ue_end_q, de_start_q, ue_end_t, de_start_t, me_start_t, me_end_t;
    int num_me_found_in_anno = 0;
    int min_exon_size = opt->min_exon_size;
    int64_t k;
    AlnParam ap = aln_param_bwa;
    int local_max_diff = (int) opt->max_diff;
    if(local_max_diff >2 ) local_max_diff = 2;
    ap.band_width = local_max_diff >1 ? local_max_diff: 1;

	uint32_t *ue_adjust_diff =new uint32_t[ n_backsearch +1 ];
	int diff_track = 0;
	for ( int ki =0; ki< n_backsearch + 1; ki++) {
		int64_t kt = ue->end_t -ki;
		int64_t kq = ue->end_q - ki;
		ubyte_t t = get_pacseq_base (pacseq, kt);
		ubyte_t q = query_seq [kq];
		if (q != t) diff_track++;
		ue_adjust_diff[ki] = diff_track;
	}
	uint32_t *de_adjust_diff = new uint32_t [  n_backsearch +1 ];
	diff_track = 0;
	for ( int ki =0; ki< n_backsearch + 1; ki++) {
		int64_t kt = de->start_t + ki;
		int64_t kq = de->start_q + ki;
		ubyte_t t = get_pacseq_base (pacseq, kt);
		ubyte_t q = query_seq [kq];
		if (q != t) diff_track++;
		de_adjust_diff[ki] = diff_track;
	}

    //TODO: need to check some place for +1 or -1
    list< pair<int64_t, int64_t> > us_sites;
    list< pair<int64_t, int64_t> > ds_sites;
	
    for (ue_end_t = min_ue_end_t; ue_end_t <= ue->end_t; ue_end_t++ ){
	//skip if the upstream pos is too close to the read start
	int breakpoint_pos_q_tmp = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
	if(breakpoint_pos_q_tmp +1 < opt->known_junc_min_anchor ) continue;
	for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
	    
	    if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != ue->strand ) ) continue;
	    if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == ue->strand ) ) continue;
			       
	    uint32_t us_type = sense_strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;
	    splice_site_t* known_ss_us= retrieve_splice_site (opt->splice_site_map, ue_end_t +1, sense_strand, us_type);
	    if(!known_ss_us) continue;
	    //check if there is a partner between ue and de
	    for(int partner_i = 0; partner_i < known_ss_us->n_partners; partner_i++ ){
		splice_site_t* partner_ss = known_ss_us->partners[partner_i];
		if(partner_ss->pos < (unsigned)de->start_t - opt->min_intron_size - min_exon_size){
		    // TODO: de->start_t is a little loose here
		    // put the pair into a list 
		    pair <int64_t, int64_t>  us_site (ue_end_t, partner_ss->pos);
		    us_sites.push_back(us_site);
		}
		
	    }
	}
    }
    
    //do the same for the downstream
    for(de_start_t = de->start_t + n_backsearch; de_start_t >= de->start_t; de_start_t--){
	//skip when the pos is too close to the end of the read
	int breakpoint_pos_q_tmp = de->start_q + de_start_t -de->start_t;
	if ( (int (seq_len) - breakpoint_pos_q_tmp) < opt->known_junc_min_anchor ) continue;
	for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {

	    if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != de->strand ) ) continue;
	    if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == de->strand ) ) continue;
			
	    uint32_t ds_type = sense_strand ? SPLICE_DONOR: SPLICE_ACCEPTOR;
	    splice_site_t* known_ss_ds= retrieve_splice_site (opt->splice_site_map, de_start_t - 1, sense_strand, ds_type);
	    if(!known_ss_ds) continue;
	    for(int partner_i = 0; partner_i < known_ss_ds->n_partners; partner_i++ ){
		splice_site_t* partner_ss = known_ss_ds->partners[partner_i];
		if(partner_ss->pos >(unsigned) min_ue_end_t + opt->min_intron_size + min_exon_size){
		    pair <int64_t, int64_t>  ds_site (partner_ss->pos, de_start_t);
		    ds_sites.push_back(ds_site);
		}
	    }
	    	    
	}
    }

    list< pair<int64_t, int64_t> >::iterator us_iter, ds_iter;
    for(us_iter = us_sites.begin(); us_iter != us_sites.end(); us_iter++ ){
	for(ds_iter = ds_sites.begin(); ds_iter != ds_sites.end(); ds_iter++ ){
	    //check the distance
	    if(us_iter->second + min_exon_size > ds_iter->first) continue;
	    for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
		
		if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != ue->strand ) ) continue;
		if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == ue->strand ) ) continue;
			    
		//check if they are on the same strand
		uint32_t us_type = sense_strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;
		splice_site_t* known_ss_us= retrieve_splice_site (opt->splice_site_map, us_iter->first +1, sense_strand, us_type);
		uint32_t ds_type = sense_strand ? SPLICE_DONOR: SPLICE_ACCEPTOR;
		splice_site_t* known_ss_ds= retrieve_splice_site (opt->splice_site_map, ds_iter->second - 1, sense_strand, ds_type);
		if(known_ss_us && known_ss_ds){
		    ue_end_t = us_iter->first;
		    ue_end_q = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
		    me_start_t = us_iter->second +1;
		    me_end_t = ds_iter->first -1 ;
		    de_start_t = ds_iter->second;
		    //de_start_q = de->start_q + n_backsearch + (de_start_t - de->start_t - n_backsearch);
		    de_start_q = de->start_q + de_start_t -de->start_t;
		    //int gap_len = de_start_q - 1 - ue_end_q -1 +1; 
		    int gap_len = de_start_q  - ue_end_q -1;
		    if(gap_len <min_exon_size) continue;
		    if(me_end_t - me_start_t +1 - gap_len > local_max_diff) continue;
		    ubyte_t *gap_seq = seq + ue_end_q + 1;
		    ubyte_t *ref_seq = (ubyte_t *) malloc(gap_len * sizeof(ubyte_t));
		    int i = 0;
		    for (k = me_start_t; i<gap_len; k++, i++)
			ref_seq[i] = get_pacseq_base (pacseq,k);
		    path_t *path = (path_t*)calloc(gap_len * 2, sizeof(path_t));
		    int path_len;
		    aln_global_core(ref_seq, gap_len, gap_seq, gap_len, &ap, path, &path_len);
		    uint32_t n_mm, n_gapo_t, n_gapo_q, n_gape_t, n_gape_q;
		    jigsaw_cal_diff (path, path_len, ref_seq, gap_len, gap_seq, gap_len,
		    &n_mm, &n_gapo_t, &n_gapo_q, &n_gape_t, &n_gape_q);
		    free (path);
		    free (ref_seq);
		    uint32_t u_an_mm = 0; 
		    uint32_t d_an_mm = 0;
		    if (ue_end_t - ue->end_t <=0) u_an_mm =  ue_adjust_diff[ ue->end_t-ue_end_t ];
		    if (de_start_t  - de->start_t >=0) d_an_mm = de_adjust_diff[de_start_t - de->start_t];
		    int diff = n_mm + n_gapo_t + n_gapo_q + n_gape_t + n_gape_q - u_an_mm - d_an_mm;
		    if (diff <= local_max_diff) {
			ue->is_last = 0;
			//create a middle exon
			jigsaw_exon_t *me= (jigsaw_exon_t *)calloc(1, sizeof(jigsaw_exon_t));
			me->n_mm= n_mm;
			me->n_gapo_t = n_gapo_t;
			me->n_gape_t = n_gape_t;
			me->n_gapo_q = n_gapo_q;
			me->n_gape_q = n_gape_q;
			me->colinear = 1;
			me->coverage = 1;
			me->strand = ue->strand;
			me->start_t= me_start_t;	
			me->end_t = me_end_t;
			me->start_q = ue_end_q + 1;
			me->end_q = de_start_q -1;
			me->hits = 0;
			me->is_first = 0;
			me->is_last = 0;
			exons->push_back(me);
			
			de->is_first =0;
			jigsaw_junction_t *us_junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
			us_junction->uexon = ue;
			us_junction->dexon = me;
			us_junction->sense_strand = sense_strand;
			us_junction->n_mm = us_junction->n_gapo_t = us_junction->n_gapo_q = 0;
			us_junction->an_mm = u_an_mm;
			us_junction->n_gape_t = us_junction->n_gape_q = 0;
			us_junction->start_q = ue_end_q;
			us_junction->end_q = us_junction->start_q +1;
			us_junction->start_t =  ue_end_t;
			us_junction->end_t = me->start_t;
			//us_junction->logistic_prob = get_splice_score(pacseq, ue_end_t, me->start_t);
			us_junction->logistic_prob = 1;
			junctions->push_back(us_junction);
			jigsaw_junction_t *ds_junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
			ds_junction->uexon = me;
			ds_junction->dexon = de;
			ds_junction->sense_strand = sense_strand;
			ds_junction->n_mm = ds_junction->n_gapo_t = ds_junction->n_gapo_q = 0;
			ds_junction->an_mm = d_an_mm;
			ds_junction->n_gape_t = ds_junction->n_gape_q = 0;
			ds_junction->start_q = me->end_q;
			ds_junction->end_q = ds_junction->start_q +1;
			ds_junction->start_t = me->end_t;
			ds_junction->end_t = de_start_t;
			//ds_junction->logistic_prob = get_splice_score(pacseq, me->end_t, de_start_t);
			ds_junction->logistic_prob = 1;
			junctions->push_back(ds_junction);
			
			num_me_found_in_anno++; 
		    }

		    
		}
	    }
		
	}
	
    }
    delete[] ue_adjust_diff;
    delete[] de_adjust_diff;
    return num_me_found_in_anno;
}

//this is to find small exons between these two anchors.
//TODO: modify this to output multi junctions later
int jigsaw_locate_junc_two_anchors_inner_exon_denovo(jigsaw_exon_t *ue, jigsaw_exon_t *de, int64_t min_ue_end_t, int64_t max_ue_end_t, int n_backsearch, uint32_t seq_len, ubyte_t *seq, const ubyte_t *pacseq,const ubyte_t *ntpac, int64_t l_pac, bwt_t *const bwt[2], const int *g_log_n, const gap_opt_t *opt, list<jigsaw_junction_t*> *junctions, list<jigsaw_exon_t*> *exons)
{
    int num_me_found_denovo = 0;
    ubyte_t *query_seq = seq;
        uint32_t *ue_adjust_diff =new uint32_t[ n_backsearch +1 ];
        int diff_track = 0;
        for ( int ki =0; ki< n_backsearch + 1; ki++) {
                int64_t kt = ue->end_t -ki;
                int64_t kq = ue->end_q - ki;
                ubyte_t t = get_pacseq_base (pacseq, kt);
                ubyte_t q = query_seq [kq];
                if (q != t) diff_track++;
                ue_adjust_diff[ki] = diff_track;
        }
        uint32_t *de_adjust_diff = new uint32_t [  n_backsearch +1 ];
        diff_track = 0;
        for ( int ki =0; ki< n_backsearch + 1; ki++) {
                int64_t kt = de->start_t + ki;
                int64_t kq = de->start_q + ki;
                ubyte_t t = get_pacseq_base (pacseq, kt);
                ubyte_t q = query_seq [kq];
                if (q != t) diff_track++;
                de_adjust_diff[ki] = diff_track;
        }

    int64_t ue_end_q, de_start_q, ue_end_t, de_start_t;
    for (ue_end_t = min_ue_end_t, ue_end_q = ue->end_q - n_backsearch; ue_end_t <= ue->end_t; ue_end_t++, ue_end_q++ )	{
	//skip if the upstream pos is too close to the read start
	int breakpoint_pos_q_tmp = ue->end_q - n_backsearch + ue_end_t - min_ue_end_t;
	if(breakpoint_pos_q_tmp +1 < opt->min_anchor ) continue;
	for (uint32_t sense_strand = 0; sense_strand < 2; ++sense_strand) {
	    
	    if( (! (opt->strand_mode & STRAND_MODE_REVERSE) ) && ( sense_strand != ue->strand ) ) continue;
	    if( (! (opt->strand_mode & STRAND_MODE_FORWARD) ) && ( sense_strand == ue->strand ) ) continue;
			
	    //always going downstream, so 
	    //upstream type:
	    uint32_t us_type = sense_strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;
	    
	    //check if this pos is a splicesite
	    bool splice_site_flag = 0;
	    
	    splice_site_t* known_ss_us= retrieve_splice_site (opt->splice_site_map, ue_end_t +1, sense_strand, us_type);
	    //note, retrieve_splice_site and jigsaw_is_splice_site_pos are tricky...pay attention to the pos used
	    if (known_ss_us) {
		    splice_site_flag =1;
	    }
	    else if (jigsaw_is_splice_site_pos ( ue_end_t +1, us_type, sense_strand, l_pac, pacseq)) {
		    splice_site_flag =1;
	    }
	    
	    if(splice_site_flag) {

		//search on de for splicesites
		uint32_t ds_type = sense_strand ? SPLICE_DONOR : SPLICE_ACCEPTOR;
		for(de_start_t = de->start_t + n_backsearch, de_start_q = de->start_q + n_backsearch; de_start_t >= de->start_t; de_start_t--, de_start_q--)	{
		     if(de_start_t - ue_end_t < 2*opt->min_intron_size + opt->min_exon_size ) continue;
			//skip when the pos is too close to the end of the read
		     breakpoint_pos_q_tmp = de->start_q + de_start_t -de->start_t;
		     if ( (int (seq_len) - breakpoint_pos_q_tmp) < opt->min_anchor ) continue;
		     int64_t length = de_start_q - ue_end_q -1;
		     if (length < opt->min_exon_size) continue;
			splice_site_flag = 0;
			splice_site_t* known_ss_ds= retrieve_splice_site (opt->splice_site_map, de_start_t - 1, sense_strand, ds_type);
			//note, retrieve_splice_site and jigsaw_is_splice_site_pos are tricky...pay attention to the pos used, use intron end here
			if (known_ss_ds) {
				splice_site_flag = 1;
			}
			else if (jigsaw_is_splice_site_pos (de_start_t - 1, ds_type, sense_strand, l_pac, pacseq)) {
				splice_site_flag = 1;
			}
			//cool to go 
			if(splice_site_flag)
			{
			    jigsaw_anchor_seq_t *anchor_seq = (jigsaw_anchor_seq_t*) calloc (1, sizeof(jigsaw_anchor_seq_t) );
			    //get splice sites
			    ubyte_t *us_ss = (ubyte_t*) calloc (2, sizeof(ubyte_t) );
			    ubyte_t *ds_ss = (ubyte_t*) calloc (2, sizeof(ubyte_t) );
			    jigsaw_get_partner_splice_site(us_type, sense_strand, us_ss);
			    jigsaw_get_partner_splice_site(ds_type, sense_strand, ds_ss);
			    anchor_seq->len = anchor_seq->full_len = length +2 +2;
			    anchor_seq->seq = (ubyte_t*)calloc(anchor_seq->len, 1);
			    
			    int j=0;
			    
			    for (int i = 0; i<2; i++,j++)
			    {
				    anchor_seq->seq[j] = us_ss[i];
			    }
			    free(us_ss);

			    for (int i = ue_end_q+1; i< de_start_q; i++,j++ )
			    {
				    anchor_seq->seq[j] = seq [i];
			    }

			    for (int i = 0; i<2; i++,j++)
			    {
				    anchor_seq->seq[j] = ds_ss[i];
			    }
			    free(ds_ss);
			    
			    anchor_seq->rseq = (ubyte_t*)calloc( anchor_seq->len, 1);
			    memcpy(anchor_seq->rseq, anchor_seq->seq, anchor_seq->len);
			    seq_reverse(anchor_seq->len, anchor_seq->rseq, opt->mode & BWA_MODE_COMPREAD);
			    anchor_seq->wid = 0;
			    anchor_seq->offset = 0;
			    anchor_seq->tid = -1;
			    anchor_seq->qual = 0;
			    //TODO: should copy the qual later
			    
			    list<jigsaw_word_hit_t*> hits;
			    jigsaw_collect_anchor_hits (bwt, anchor_seq, g_log_n,
				    l_pac, pacseq, ntpac, opt->max_word_occ, opt, &hits);
			    free (anchor_seq->seq);
			    free (anchor_seq->rseq);
			    free (anchor_seq);
			    if(hits.size() == 0 ) continue;
			    jigsaw_sort_word_hits_by_strand_pos_t(&hits);
			    int64_t left_bound_t = ue_end_t + opt->min_intron_size;
			    int64_t right_bound_t = de_start_t - opt->min_intron_size - length;
			    
			    list<jigsaw_word_hit_t*>::iterator iter = hits.begin();
			    
			    for ( ; iter != hits.end(); iter++){
			
				jigsaw_word_hit_t *piter = *iter;
				int64_t p_pos_t = piter->pos_t;
				if (piter->strand == 0) {
				    if ( p_pos_t <right_bound_t ) {
					if (p_pos_t > left_bound_t ){
					    ue->is_last = 0;
					    //create a middle exon
					    jigsaw_exon_t *me= (jigsaw_exon_t *)calloc(1, sizeof(jigsaw_exon_t));
					    me->n_mm= me->n_gapo_t = me->n_gape_t = me->n_gapo_q = me->n_gape_q = 0;
					    me->colinear = 1; me->coverage = 1;
					    me->strand = ue->strand;
					    me->start_t = piter->pos_t +2;
					    me->end_t =  me->start_t + length -1;
					    me->start_q = ue_end_q+1;
					    me->end_q = ue_end_q + length;
					    me->hits = 0;
					    me->is_first = 0;
					    me->is_last = 0;

					    //modify the downstream exon, no need
					    //de->start_q = de_start_q;
					    //de->start_t = de_start_t;
					    de->is_first =0;
					    // the first junction

					    jigsaw_junction_t *us_junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
					    us_junction->uexon = ue;
					    us_junction->dexon = me;
					    us_junction->sense_strand = sense_strand;
					    us_junction->n_mm = us_junction->n_gapo_t = us_junction->n_gapo_q = 0;
					    us_junction->an_mm = 0;
					    if (ue_end_t - ue->end_t <=0) us_junction->an_mm = ue_adjust_diff[ ue->end_t-ue_end_t ];
					    us_junction->n_gape_t = us_junction->n_gape_q = 0;
					    us_junction->start_q = ue_end_q;
					    us_junction->end_q = us_junction->start_q +1;
					    us_junction->start_t = ue_end_t;
					    us_junction->end_t = me->start_t;
				    us_junction->logistic_prob = get_splice_score(pacseq, ue_end_t, me->start_t);

					    junctions->push_back(us_junction);
					    
					    //the 2nd exon
					    jigsaw_junction_t *ds_junction = (jigsaw_junction_t*) calloc (1, sizeof(jigsaw_junction_t));
					    ds_junction->uexon = me;
					    ds_junction->dexon = de;
					    ds_junction->sense_strand = sense_strand;
					    ds_junction->n_mm = ds_junction->n_gapo_t = ds_junction->n_gapo_q = 0;
					    ds_junction->an_mm = 0;
					    if (de_start_t  - de->start_t >=0) ds_junction->an_mm = de_adjust_diff[de_start_t - de->start_t];
					    ds_junction->n_gape_t = ds_junction->n_gape_q = 0;
					    ds_junction->start_q = me->end_q;
					    ds_junction->end_q = ds_junction->start_q +1;
					    ds_junction->start_t = me->end_t;
					    ds_junction->end_t = de_start_t;
				    ds_junction->logistic_prob = get_splice_score(pacseq, me->end_t, de_start_t);

					    junctions->push_back(ds_junction);
					    exons->push_back(me);
					    num_me_found_denovo++;
					}
				    }
				    else break;
				}
				else break;
			    }
		    
			    for (iter = hits.begin() ; iter != hits.end(); iter++)
			    {
				    jigsaw_word_hit_t *piter = *iter;
				    free(piter);
				    
			    }

		    }
			
		}
	}
    }
	    
    }
    delete[] ue_adjust_diff;
    delete[] de_adjust_diff;
    return num_me_found_denovo;

}


/*
 * TODO: 1. it does not allow indels during local search of splice sites now, need to be extended later
 */
void jigsaw_locate_junc_two_anchors (jigsaw_exon_t *ue, jigsaw_exon_t *de,
                ubyte_t *seq, uint32_t seq_len, const ubyte_t *pacseq, const ubyte_t *ntpac, int64_t l_pac, bwt_t *const bwt[2], const int *g_log_n, const gap_opt_t *opt,
                list<jigsaw_junction_t*> *junctions, list<jigsaw_exon_t*> *exons)

{
	//TODO: there might be over estimate when we need to come back to find the splice sites, and there is mismatch in the region
	//specified by max_overhang, should be relatively rare though

	int n_backsearch = opt->max_overhang;

	//when sequences extends beyond splice sites by chance, we can increase the overhang
	if (n_backsearch < ue->end_q - de->start_q + 1) n_backsearch = ue->end_q - de->start_q + 1;

	if (ue->end_q - ue->start_q + 1 - n_backsearch <= 4)
		n_backsearch = ue->end_q - ue->start_q + 1 - 4;

	int gap_start_q = ue->end_q - n_backsearch + 1;
	int gap_len = de->start_q - ue->end_q - 1 + 2 * n_backsearch;
	if (gap_len <= 0) return;

	if (gap_len >= (int) seq_len - gap_start_q) gap_len = (int) seq_len - gap_start_q;

	if (de->start_t + 2*n_backsearch - gap_len -ue->end_t <= 0 ) return;
	//in case gap in the target sequence is too big, this may lead to 
	//problem of "negative intron length", to ensure de_start_t > ue_end_t in the loop below
	//the ineq. above was derived to bypass this occassion. 
	//bug found 0913, 
	
	//TODO: the start and end of gap needs to be checked more carefully

	AlnParam ap = aln_param_bwa;
	//TODO: need more sophisticated rules; the bound is too loose here
	ap.band_width =  int (opt->max_diff );
	if (ap.band_width <1) ap.band_width = 1;
	if( (!opt->splice_site_map) && ap.band_width>2) ap.band_width = 2;

	int64_t min_ue_end_t = ue->end_t - n_backsearch;
	int64_t max_ue_end_t = min_ue_end_t + gap_len - 1;

	//ubyte_t *gap_seq = seq + ue->end_q - n_backsearch + 1;
	ubyte_t *ref_seq = (ubyte_t *) malloc(gap_len * sizeof(ubyte_t));

	int num_junc_found_denovo = 0;
	int num_junc_found_in_anno = 0;
	
	if(opt->splice_site_map)// if junction database is provided, check it first
	{
	    num_junc_found_in_anno = jigsaw_locate_junc_two_anchors_with_anno(ue, de, min_ue_end_t, max_ue_end_t, n_backsearch, seq_len, gap_len, seq, ref_seq, pacseq, &ap, opt, junctions);
	    //if(opt->non_denovo_search) {free(ref_seq); return;}
	    //stop here if both splice_site_map and non_denovo_search

	}
	if((opt->non_denovo_search == 0) && de->start_t - 1 - ue->end_t < opt->max_intron_size) 

        {
	    num_junc_found_denovo = jigsaw_locate_junc_two_anchors_denovo(ue, de, min_ue_end_t, max_ue_end_t, n_backsearch, seq_len, gap_len, seq, ref_seq, pacseq, &ap, opt, junctions);
	
       }
       free (ref_seq);
       
     // only do this if no junction was found above 
     // or there is a junction db
     if(opt->splice_site_map ||( num_junc_found_denovo ==0 && num_junc_found_in_anno == 0) ) {
	int num_me_found_in_anno = 0, num_me_found_denovo = 0;	
	if (de->start_t - 1 - ue->end_t > 2 * opt->max_intron_size + opt->min_anchor) return;
	 if(opt->splice_site_map) {
	     num_me_found_in_anno = jigsaw_locate_junc_two_anchors_inner_exon_with_anno(ue, de, min_ue_end_t, max_ue_end_t, n_backsearch, seq_len, seq, pacseq, ntpac, l_pac,  bwt, g_log_n, opt, junctions, exons);
	 }
	 if (opt->non_denovo_search == 0){
	     num_me_found_denovo = jigsaw_locate_junc_two_anchors_inner_exon_denovo(ue, de, min_ue_end_t, max_ue_end_t, n_backsearch, seq_len, seq, pacseq, ntpac, l_pac,  bwt, g_log_n, opt, junctions, exons);
	 }
      }
}



void jigsaw_pair_exons (bwa_seq_t *seq, list<jigsaw_exon_t*> *exons, const gap_opt_t *opt,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac, bwt_t *const bwt[2], const int *g_log_n, list<jigsaw_junction_t*> *junctions)
{
	int word_size = seq->word_size;
	jigsaw_exon_t *ue, *de, *e;
	list<jigsaw_exon_t*>::iterator iter = exons->begin();
	//uint32_t strand = (*iter)->strand;
	int64_t gap;

	//fprintf (stderr, "pair exons\n");

	//fprintf (stderr, "locate junction with two anchors\n");

	//locate junctions using two anchors
	int diff, wid_diff;
	int n = exons->size(), i, j;
	for (iter = exons->begin(), i=0; i<n; ++iter, ++i) {
		ue = *iter;
		diff = ue->n_mm + ue->n_gapo_t + ue->n_gape_t + ue->n_gapo_q + ue->n_gape_q;
		if (ue->colinear == 0 || diff > (int)opt->max_diff) continue;
		list<jigsaw_exon_t*>::iterator iter2 = iter;
		j = i + 1;
		for (iter2++; j<n; ++iter2, ++j) {
			de = *iter2;
			if(ue->strand != de->strand ) continue;
			diff = de->n_mm + de->n_gapo_t + de->n_gape_t + de->n_gapo_q + de->n_gape_q;
			wid_diff = de->hits->front ()->wid - ue->hits->back()->wid;
			gap = de->start_q - ue->end_q - 1;

			//note that some words might be repetitive so not used as seeds, but they might be recovered by extension
			//of neighboring seeds, so gap might be smaller
			if (de->colinear == 0 || diff > (int)opt->max_diff || wid_diff <= 0
					 || de->start_t <= ue->end_t) continue;
			if ( (!opt->splice_site_map) && gap > 2 * word_size) continue;

			if ( (!opt->splice_site_map) && de->start_t - 1 - ue->end_t > 2 * opt->max_intron_size) break;
			//two times of max_intron size, since we may have an exon between them
			jigsaw_locate_junc_two_anchors (ue, de, ue->strand ? seq->rseq : seq->seq, seq->len, pacseq, ntpac, l_pac, bwt, g_log_n, opt, junctions, exons);


	
		}
	}

	//fprintf (stderr, "locate junction with one anchors\n");


	if (opt->splice_site_map || opt->single_anchor_search) {
		int n = exons->size(), i;
		for (i = 0, iter = exons->begin(); i < n; ++iter, ++i) {
			e = *iter;
			if (e->is_first) {

				//locate upstream junction
				uint32_t direction = 1;

				//TODO: some times the other end are long but has sequencing errors
				//those will be missed at the moment
				if( (e->start_q + opt->max_overhang >= opt->min_anchor && e->start_q < 2 * word_size) || opt->splice_site_map )
				{

					//fprintf (stderr, "extend exon %d start at %d on the left\n", e->exon_id, e->start_q);
					jigsaw_locate_junc_one_anchor (e, seq, direction,
							opt, l_pac, pacseq, ntpac, bwt, g_log_n, junctions, exons);

				}
			}

			if (e->is_last) {


				//locate downstream junction
				uint32_t direction = 0;

				if ( (seq->len - e->end_q - 1 + opt->max_overhang >= opt->min_anchor && seq->len - e->end_q - 1 < 2* word_size) || opt->splice_site_map)
					//allow more single anchor searching than before.(overhang)
				{
					//fprintf (stderr, "extend exon %d end at %d on the right\n", e->exon_id, e->end_q);

					jigsaw_locate_junc_one_anchor (e, seq, direction,
							opt, l_pac, pacseq, ntpac, bwt, g_log_n, junctions, exons);

					}
			}
		}//loop exon
	}
}


/* compare hits by strand, and then by pos_t
 */
bool jigsaw_junction_comp_strand_pos_t (const jigsaw_junction_t *a, const jigsaw_junction_t *b)
{
	if (a->uexon->strand < b->uexon->strand) return true;
	else if (a->uexon->strand > b->uexon->strand)return false; //positive strand (0) first

	if (a->sense_strand < b->sense_strand) return true;
	else if (a->sense_strand > b->sense_strand) return false;

	//the same strand compare position
	int64_t sa = a->uexon->start_t;
	int64_t sb = b->uexon->start_t;
	if (sa < sb ) return true;
	else if (sa > sb ) return false;
	

	if (a->start_t <  b->start_t) return true;
	else if (a->start_t >  b->start_t) return false;
	else return (a->end_t < b->end_t);
}

/* sort hits (ascendingly) by strand and position on the genome
 */
void jigsaw_sort_junctions_by_strand_pos_t (list<jigsaw_junction_t*> *junctions)
{
	junctions->sort (jigsaw_junction_comp_strand_pos_t);
}

void jigsaw_uniq_junctions ( list<jigsaw_junction_t*> *junctions )
{
   list<jigsaw_junction_t*>::iterator iter = junctions->begin();
   list<jigsaw_junction_t*>::iterator last_iter ;
   iter++;
   while(iter != junctions->end()){
       last_iter = iter;
       last_iter--;
       jigsaw_junction_t *p = *iter;
       jigsaw_junction_t *last_p = *last_iter;
       //if(p->strand == last_p->strand && p->uexon == last_p->uexon && p->dexon == last_p->dexon && p->start_t == last_p->start_t && p->end_t == last_p->end_t) {
	// TODO : need to think about the strand again
       if(p->uexon == last_p->uexon && p->dexon == last_p->dexon && p->start_t == last_p->start_t && p->end_t == last_p->end_t) {
	   free(p);
	   iter = junctions->erase(iter);
       }
       else{
	   iter++; 
       }
   }
}

/* 
 * alignment will be saved in aln
 */
void jigsaw_concat_junctions (list <jigsaw_junction_t*> *junctions, int len,
		uint32_t global, list<jigsaw_spliced_aln_t*> *aln)
{
	//int i, j;
	//int _n_aln = *n_aln, _n_aln_old = *n_aln, _max_aln = *max_aln;

	if (junctions->size() == 0) return;
	jigsaw_sort_junctions_by_strand_pos_t (junctions);
	//sometimes there are duplicate junctions from both with_anno and denovo search
	jigsaw_uniq_junctions (junctions);

	list <jigsaw_spliced_aln_t*>::iterator aln_iter;
	//get the first junction
	jigsaw_junction_t *p = junctions->front();
	jigsaw_junction_t *curr_p; //this is going to be tmp pointer to a jigsaw_junction_t
	jigsaw_junction_t *last_p; //I am another tmp pointer

	//create the first alignment and add the "first" junction
	jigsaw_spliced_aln_t *curr_aln = (jigsaw_spliced_aln_t *) calloc (1, sizeof (jigsaw_spliced_aln_t));
	curr_aln->junctions = new list<jigsaw_junction_t*>;
	curr_aln->junctions->push_back (p);
	aln->push_back (curr_aln);

	//as an initiation, put the junctions starting at the same uexon as the "first" junction into the list<aln>
	list <jigsaw_junction_t*>::iterator iter = junctions->begin(); ++iter;
	//loop start from the second junction
	for (; iter != junctions->end(); iter++) {
	    curr_p = *iter;
	    if ( p->uexon == curr_p->uexon && curr_p->sense_strand == p->sense_strand ){
		//check if this junction (curr_p) has the same start exon as the "first" junction (p)
		//make a new branch start
		curr_aln = (jigsaw_spliced_aln_t *) calloc (1, sizeof (jigsaw_spliced_aln_t));
		curr_aln->junctions = new list<jigsaw_junction_t*>;
		curr_aln->junctions->push_back (curr_p);
		aln->push_back (curr_aln);
	    }
	    else {
		//end the initiation
		break;
	    }
	}
		
	//record the position of the first alignment
	//--first_aln_iter;
	//continue from the current iter of junctions
	int previous_aln_size = int(aln->size() );
	//record the size of aln, only search the first n alignment in each loop, these are the
	//alignment before the new branches
	bool *destroy_flag =  new bool[previous_aln_size];
	// whether the current aln should be destroy 1: yes, 0 no
	for(int i=0; i<previous_aln_size; i++){
	    destroy_flag[i] = 0;
	}
	list <jigsaw_junction_t*>::iterator last_iter = iter; last_iter--;	
	bool no_need_to_update_aln_list = 1;
	while (iter != junctions->end()) {
		curr_p = *iter;
		//TODO: could do some speedup since it's inside one cluster, now the code is for general use
		//could stop or do some recovery if the aln is not conitnuous. but list<aln> comes from all the clusters, so its big and everytime we go from the first aln, might do some optimization here. 
		//leave it like this for now Apr 11th WJ
		
		//since junctions are ordered according to the start position
		//we should seek the possibility of extension in the existing alignment
		//we also need to make branch if the junction "overlaps" with the last one, in the multiple output mode
		// extend if uexon of the current junc is the last exon in aln, 
		// first make a copy of this alignment, push it into list<aln>, 
		// then this junc will be pushed into aln->junc,
		// there will be only two cases 
		last_p = *last_iter;
		if(last_p->uexon != curr_p->uexon && (!no_need_to_update_aln_list) ){
		    //update list<aln> since the junction is moving forward
		    //for now, no_need_to_update_aln_list will only happen in the first loop
		    no_need_to_update_aln_list = 1;
		    aln_iter = aln->begin();
		    for (int i = 0; i<previous_aln_size;  i++){
			if(destroy_flag[i]){
			    //destroy this aln in list<aln>
			    jigsaw_spliced_aln_t *aln_p_tmp=*aln_iter;
			    delete(aln_p_tmp->junctions);
			    free (aln_p_tmp);
			    aln_iter = aln->erase(aln_iter);
			    //refer to list::erase, aln_iter is the iterator to the next pos already
			}
			else{
			    aln_iter++;
			}
		    }
		    previous_aln_size = int(aln->size() );
		    //update the previous_aln_size
		    delete[] destroy_flag;
		    destroy_flag =  new bool[previous_aln_size];
		    //update the destroy_flag
		    for(int i=0; i<previous_aln_size; i++){
			destroy_flag[i] = 0;
		    }
		}
		
		bool extended = 0;//for curr junction
		aln_iter = aln->begin();
		for (int i = 0; i < previous_aln_size; aln_iter++, i++) {
		    //the alns after previous_aln_size are newly created copies
			curr_aln = *aln_iter;
			//get the last junction in the current alignment
			jigsaw_junction_t *last = curr_aln->junctions->back();
			if (last->dexon == curr_p->uexon && last->sense_strand == curr_p->sense_strand && last->start_q < curr_p->start_q) { 
			    //the uexon of this junction (curr_p) matches the upstream dexon in aln
			    //make a copy of curr_aln
			    jigsaw_spliced_aln_t *curr_aln_copy = (jigsaw_spliced_aln_t *) calloc (1, sizeof (jigsaw_spliced_aln_t));
			    curr_aln_copy->junctions = new list<jigsaw_junction_t*>;
			    list <jigsaw_junction_t*>::iterator iter_copy = curr_aln->junctions->begin();
			    //copy the junctions to curr_aln_copy
			    for (; iter_copy != curr_aln->junctions->end(); iter_copy++) {
				jigsaw_junction_t *p_copy = *iter_copy;
				curr_aln_copy->junctions->push_back(p_copy);
			    }
			    //extend: push the current junction curr_p
			    curr_aln_copy->junctions->push_back(curr_p);
			    //push this copy into the end of list<aln>
			    aln->push_back(curr_aln_copy);
			    no_need_to_update_aln_list = 0;
			    if ( last->dexon->end_q != len-1 ) destroy_flag [i] = 1;// this aln will be destroyed later
				//if the last exon reached the boundary, keep it as a candidate aln
			    extended = 1;
			    
			}
		}

		if( extended ==0 || curr_p->uexon->start_q == 0 ){
		    // it can not be extended, 
		    // or the extension step reached the boundary, this is to fix the bug discovered by Rahul, shown below, when there can be two solutions, the previous concate function only reported the first solution, since junction 2, which although could form a single aln, is connected into the longer aln. 
		    // =====    1       ====     2       ======
		    //            ==========     2       ======
		    //create a new aln
		    curr_aln = (jigsaw_spliced_aln_t *) calloc (1, sizeof (jigsaw_spliced_aln_t));
		    curr_aln->junctions = new list<jigsaw_junction_t*>;
		    curr_aln->junctions->push_back(curr_p);
		    aln->push_back(curr_aln);
		    no_need_to_update_aln_list = 0;
		}
		
		//move to the next junction
		last_iter++;
		iter++;

	}
	
	if(!no_need_to_update_aln_list) {
	//update list<aln> 
	    aln_iter = aln->begin();
	    for (int i = 0; i<previous_aln_size;  i++){
		if(destroy_flag[i]){
		    //destroy this aln in list<aln>
		    jigsaw_spliced_aln_t *aln_p_tmp=*aln_iter;
		    delete(aln_p_tmp->junctions);
		    free (aln_p_tmp);
		    aln_iter = aln->erase(aln_iter);
		    //refer to list::erase, aln_iter is the iterator to the next pos already
		}
		else{
		    aln_iter++;
		}
	    }
	}

	delete[] destroy_flag;

	//evaluate alignments
	//destroy duplicate alns in the mean time
	list <jigsaw_spliced_aln_t*>::iterator last_aln_iter;
	for (aln_iter = aln->begin(); aln_iter != aln->end (); ++aln_iter) {
		curr_aln = *aln_iter;

		//check if this aln is a duplicate one, this could come from inner exon search
		if(aln_iter != aln->begin()){
		    bool is_identical = 1;
		    last_aln_iter = aln_iter;
		    last_aln_iter--;
		    jigsaw_spliced_aln_t *last_aln = *last_aln_iter;
		    if(curr_aln->junctions->size() == last_aln->junctions->size()){
			last_iter = last_aln->junctions->begin();
			iter = curr_aln->junctions->begin();
			for (; iter != curr_aln->junctions->end(); ++iter, ++last_iter){
			    last_p = *last_iter;
			    curr_p = *iter;
			    if(last_p->start_t != curr_p->start_t || last_p->end_t != curr_p->end_t || last_p->uexon->strand != curr_p->uexon->strand){
				is_identical = 0;
				break;
			    }
			}
		    }
		    else {
			is_identical = 0;
		    }

		    if(is_identical){
			//destroy curr_aln
			delete(curr_aln->junctions);
			free(curr_aln);
			aln_iter = aln->erase( aln_iter );
			aln_iter--;
			continue;
		    }
		}

		jigsaw_junction_t *first = curr_aln->junctions->front();
		jigsaw_junction_t *last = curr_aln->junctions->back();

		curr_aln->strand = first->uexon->strand;
		curr_aln->sense_strand = first->sense_strand;

		curr_aln->start_t = first->uexon->start_t;
		curr_aln->start_q = first->uexon->start_q;
		curr_aln->end_t = last->dexon->end_t;
		curr_aln->end_q = last->dexon->end_q;

		curr_aln->n_mm = first->uexon->n_mm;
		curr_aln->n_gapo_t = first->uexon->n_gapo_t;
		curr_aln->n_gape_t = first->uexon->n_gape_t;
		curr_aln->n_gapo_q = first->uexon->n_gapo_q;
		curr_aln->n_gape_q = first->uexon->n_gape_q;
		curr_aln->diff = 0;
		curr_aln->logistic_prob = 0.0;

		for (iter = curr_aln->junctions->begin(); iter != curr_aln->junctions->end(); ++iter) {
			p = *iter;
			curr_aln->n_mm += p->n_mm;
			curr_aln->n_gapo_t += p->n_gapo_t;
			curr_aln->n_gape_t += p->n_gape_t;
			curr_aln->n_gapo_q += p->n_gapo_q;
			curr_aln->n_gape_q += p->n_gape_q;

			curr_aln->n_mm += p->dexon->n_mm;
			// n_mm might be overestimated since there redundant mm counts between junction and de ue
			curr_aln->n_gapo_t += p->dexon->n_gapo_t;
			curr_aln->n_gape_t += p->dexon->n_gape_t;
			curr_aln->n_gapo_q += p->dexon->n_gapo_q;
			curr_aln->n_gape_q += p->dexon->n_gape_q;

			curr_aln->logistic_prob += p->logistic_prob;

			//adjust the mismatch with an_mm in junction
			//do it at the end of the loop to avoid negative numbers
			curr_aln->n_mm -= p->an_mm;
		}

		curr_aln->diff = curr_aln->n_mm + curr_aln->n_gapo_t + curr_aln->n_gape_t
				+ curr_aln->n_gapo_q + curr_aln->n_gape_q;
		curr_aln->logistic_prob = curr_aln->logistic_prob / curr_aln->junctions->size();

		if (global) {
			curr_aln->diff += curr_aln->start_q;
			curr_aln->diff += (len - 1 - curr_aln->end_q);

			//avoid end truncation for global alignment
			curr_aln->start_t -= curr_aln->start_q;
			curr_aln->end_t += (len - 1 - curr_aln->end_q);
			curr_aln->start_q = 0; curr_aln->end_q = len - 1;
		}
	}
	//return aln;
}


/*compare spliced alignment by score (ascending)
 */

bool jigsaw_spliced_aln_comp_diff (const jigsaw_spliced_aln_t *a, const jigsaw_spliced_aln_t *b)
{
	//compare diff
	if (a->diff < b->diff) return true;
	else if (a->diff > b->diff) return false;

	//the number of diff is the same, then compare gap
	int gapo_a = a->n_gapo_t + a->n_gapo_q;
	int gapo_b = b->n_gapo_t + b->n_gapo_q;

	if (gapo_a < gapo_b) return true;
	else if (gapo_a > gapo_b) return false;

	int gape_a = a->n_gape_t + a->n_gape_q;
	int gape_b = b->n_gape_t + b->n_gape_q;

	if (gape_a < gape_b) return true;
	else if (gape_a > gape_b) return false;

	//compare intron size, favor smaller intron
	//int intron_size_a = a->end_t - a->start_t;
	//int intron_size_b = b->end_t - b->start_t;
	//if (intron_size_a < intron_size_b) return true;
	//else if (intron_size_a > intron_size_b) return false;

	//compare junction numbers, favor fewer junctions
	//return a->junctions->size() < b->junctions->size() ? true: false;
	return a->logistic_prob > b->logistic_prob ?  true: false;
}


/*compare spliced alignment by score (ascending)
 *
 */
void jigsaw_sort_spliced_aln_by_diff (list<jigsaw_spliced_aln_t*> *aln)
{
	//sort all clumps
	aln->sort (jigsaw_spliced_aln_comp_diff);
}

void jigsaw_search_exonic_aln (jigsaw_spliced_aln_cluster_t *clust, bwa_seq_t *seq, const gap_opt_t *opt, list<jigsaw_spliced_aln_t*> *aln)
{
	list <jigsaw_exon_t *>::iterator iter;
	for (iter = clust->exons->begin(); iter != clust->exons->end (); ++iter) {
              jigsaw_exon_t *p = *iter;
	      if(p->start_q == 0 && p->end_q == seq->len -1)
	      {
		//check if this pos_t already recorded
		bool aln_checker = 0;
		list<jigsaw_spliced_aln_t*>::iterator aln_iter;
		for (aln_iter = aln->begin(); aln_iter != aln->end (); aln_iter++){
			jigsaw_spliced_aln_t *check_aln = *aln_iter;
			if ( check_aln->junctions == NULL and check_aln->start_t == p->start_t ) {
				aln_checker = 1;
				break;
			}
		}
		if ( aln_checker ==0 ){
			  // make a new exonic alignment, use the same data structure, but no junctions
			  jigsaw_spliced_aln_t *curr_aln = (jigsaw_spliced_aln_t *) calloc (1, sizeof (jigsaw_spliced_aln_t));
			  // init' this aln
			  curr_aln->junctions = NULL;
			  curr_aln->strand = p->strand;
			  curr_aln->sense_strand = 2;// cannot be determined from data
			  curr_aln->start_t = p->start_t;
			  curr_aln->start_q = p->start_q;
			curr_aln->end_t = p->end_t;
			curr_aln->end_q = p->end_q;
			curr_aln->n_mm = p->n_mm;
			curr_aln->n_gapo_t = p->n_gapo_t;
			curr_aln->n_gape_t = p->n_gape_t;
			curr_aln->n_gapo_q = p->n_gapo_q;
			curr_aln->n_gape_q = p->n_gape_q;
			curr_aln->diff =curr_aln->n_mm + curr_aln->n_gapo_t + curr_aln->n_gape_t + curr_aln->n_gapo_q + curr_aln->n_gape_q;
			curr_aln->logistic_prob = 1.0;
			aln->push_back(curr_aln);			  
		}
	      }
	}
			     
    
}

void jigsaw_search_spliced_aln (jigsaw_spliced_aln_cluster_t *clust, bwa_seq_t *seq,
		const gap_opt_t *opt, uint32_t extend_exact, uint32_t global,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,bwt_t *const bwt[2], const int *g_log_n,  list<jigsaw_spliced_aln_t*> *aln)
{
	//here we assume all hits, and exons are on the same strand
	//jigsaw_spliced_aln_t *aln = 0;
	//int n_aln = 0, max_aln = 0;

	//if (jigsaw_spliced_aln_cluster_unique_hits (clust, seq->n_words) < (uint32_t) seq->n_words/2)
	//	return;

	jigsaw_sort_word_hits_by_diagonal (clust->hits);

	//fprintf (stderr, "\ngroup hits to exon\n");

	//fprintf(stderr, "seq name=%s, cluster_id=%d, %d hits\n", seq->name, clust->cluster_id, clust->hits->size());
	clust->exons = new list<jigsaw_exon_t *>;
	jigsaw_group_hits_to_exons (clust->hits, opt->max_diff, seq->word_size, seq->n_words, clust->exons);

	if (clust->exons->size()==0) return;

	//int n_exons = clust->exons->size();
	//free (clust->hits); clust->hits = 0;

	//fprintf (stderr, "extend exon\n");

	bool stop_if_neg_score = false;
	jigsaw_extend_exons (clust->exons, seq, l_pac, pacseq, ntpac, extend_exact, stop_if_neg_score, opt->max_diff);

	//fprintf (stderr, "refine exon\n");

	jigsaw_refine_exon_aln (clust->exons, seq, l_pac, pacseq, ntpac, opt->max_diff);
	jigsaw_sort_exons_by_pos_t (clust->exons);

	//fprintf (stderr, "pair exon\n");

	clust->junctions = new list<jigsaw_junction_t *>;
	jigsaw_pair_exons (seq, clust->exons, opt, l_pac, pacseq, ntpac, bwt, g_log_n, clust->junctions);
	//int n_junctions = clust->junctions->size();

	if (clust->junctions->size() > 0)
		jigsaw_concat_junctions (clust->junctions, seq->len, global, aln);

	//return aln;
}


//This is a very simplified version that specifies only the position of exon/introns
uint32_t *spliced_aln2cigar32_simple (const jigsaw_spliced_aln_t *aln, int *n_cigar)
{
	int n, len;
	uint32_t *cigar;
	//unsigned char last_type;

	jigsaw_junction_t *p, *q;
	int n_junctions = aln->junctions->size ();
	if (n_junctions == 0) {
		*n_cigar = 0;
		return 0;
	}

	*n_cigar = n_junctions *2 + 1;
	cigar = (uint32_t*)malloc(*n_cigar * sizeof (uint32_t));

	list<jigsaw_junction_t *>::iterator iter, iter0;

	for (n=0, iter = aln->junctions->begin(); iter != aln->junctions->end(); ++iter) {
		p = *iter;

		//the upstream exon
		if (iter == aln->junctions->begin()) //first exon
			len = p->start_t - p->uexon->start_t + 1;
		else {
			iter0 = iter; --iter0;
			q = *iter0;
			len = p->start_t - q->end_t + 1;
		}
		cigar[n++] = len << 4 | CIGAR_OP_M;

		//junction
		len = p->end_t - p->start_t - 1;
		cigar[n++] = len << 4 | CIGAR_OP_N;
	}

	//last exon
	--iter; p = *iter;
	len = p->dexon->end_t - p->end_t + 1;
	cigar[n++] = len << 4 | CIGAR_OP_M;

	return cigar;
}



/*This only deals with exonic alignment from jigsaw*/
uint32_t *exonic_aln2cigar32(const jigsaw_spliced_aln_t *aln, ubyte_t *seq,
                int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,
		                int band_width, int *n_cigar)
{
	int max_cigar, _n_cigar = 0, n_exon_cigar;
	uint32_t *cigar, *exon_cigar;
	AlnParam ap = aln_param_bwa;
	ap.band_width = band_width > 1? band_width: 1; //override band width
	max_cigar = *n_cigar = 1;
	cigar = (uint32_t*)malloc(*n_cigar * sizeof (uint32_t));
	int len_t, len_q, l, path_len;
	int64_t k;
	len_t = aln->end_t - aln->start_t +1;
	len_q = aln->end_q - aln->start_q +1;
	ubyte_t* ref_seq = (ubyte_t*)calloc(len_t, 1);
	for (k = aln->start_t, l = 0; k < aln->start_t + len_t && k < l_pac; ++k)
	        ref_seq[l++] = get_pacseq_base (pacseq, k);
	path_t* path = (path_t*)calloc(len_t+len_q, sizeof(path_t));
	aln_global_core(ref_seq, len_t, (ubyte_t*)(seq + aln->start_q), len_q, &ap, path, &path_len);
	exon_cigar = aln_path2cigar32(path, path_len, &n_exon_cigar);
	
	if (max_cigar <= _n_cigar + n_exon_cigar) {
	    max_cigar = (_n_cigar + n_exon_cigar) << 1;
	    cigar = (uint32_t*)realloc(cigar, sizeof(uint32_t) * max_cigar);
	}
	for (int i = 0; i < n_exon_cigar; ++i)
	    cigar[_n_cigar++] = exon_cigar[i];
	free (ref_seq); free (path); free(exon_cigar);
	*n_cigar = _n_cigar;
	return cigar;

				    
}

jigsaw_cigar_t *jigsaw_exonic_aln2cigar(const jigsaw_spliced_aln_t *aln, ubyte_t *seq,
int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,
		int band_width, int *n_cigar)
{
	uint32_t *cigar32;
	jigsaw_cigar_t *cigar;
	int i;
	cigar32 = exonic_aln2cigar32(aln, seq, l_pac, pacseq, ntpac, band_width, n_cigar);
	cigar = (jigsaw_cigar_t*)cigar32;
	for (i = 0; i < *n_cigar; ++i)
	cigar[i] = __cigar_create( (cigar32[i]&0xf), (cigar32[i]>>4) );
	return cigar;
}


/*This is the extended version that realign each exon to identify mismatches*/
uint32_t *spliced_aln2cigar32(const jigsaw_spliced_aln_t *aln, ubyte_t *seq,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,
		int band_width, int *n_cigar)
{
	int n, len, max_cigar, _n_cigar = 0, n_exon_cigar;
	uint32_t *cigar, *exon_cigar;
	//unsigned char last_type;

	AlnParam ap = aln_param_bwa;
	ap.band_width = band_width > 1 ? band_width : 1 ; //override band width

	jigsaw_junction_t *p, *q;
	int n_junctions = aln->junctions->size ();
	if (n_junctions == 0) {
		*n_cigar = 0;
		return 0;
	}

	max_cigar = *n_cigar = n_junctions *2 + 1;
	cigar = (uint32_t*)malloc(*n_cigar * sizeof (uint32_t));

	list<jigsaw_junction_t *>::iterator iter, iter0;
	int len_t, len_q, l, path_len;
	int64_t start_t, start_q, k;

	for (n=0, iter = aln->junctions->begin(); iter != aln->junctions->end(); ++iter) {
		p = *iter;

		//the upstream exon
		if (iter == aln->junctions->begin()) {//first exon
			len_t = p->start_t - aln->start_t + 1;
			len_q = p->start_q - aln->start_q + 1;
			start_t = aln->start_t;
			start_q = aln->start_q;
			//len_t = p->start_t - p->uexon->start_t + 1;
			//len_q = p->start_q - p->uexon->start_q + 1;
			//start_t = p->uexon->start_t;
			//start_q = p->uexon->start_q;
		}
		else {
			iter0 = iter; --iter0;
			q = *iter0;
			len_t = p->start_t - q->end_t + 1;
			len_q = p->start_q - q->end_q + 1;
			start_t = q->end_t; start_q = q->end_q;
		}

		//get ref seq
		ubyte_t* ref_seq = (ubyte_t*)calloc(len_t, 1);
		for (k = start_t, l = 0; k < start_t + len_t && k < l_pac; ++k)
			ref_seq[l++] = get_pacseq_base (pacseq, k);

		path_t* path = (path_t*)calloc(len_t+len_q, sizeof(path_t));

		aln_global_core(ref_seq, len_t, (ubyte_t*)(seq+start_q), len_q, &ap, path, &path_len);
		exon_cigar = aln_path2cigar32(path, path_len, &n_exon_cigar);

		if (max_cigar <= _n_cigar + n_exon_cigar) {
			max_cigar = (_n_cigar + n_exon_cigar) << 1;
			cigar = (uint32_t*)realloc(cigar, sizeof(uint32_t) * max_cigar);
		}

		for (int i = 0; i < n_exon_cigar; ++i)
			cigar[_n_cigar++] = exon_cigar[i];

		free (ref_seq); free (path); free(exon_cigar);

		//intron
		len = p->end_t - p->start_t - 1;
		xassert (len > 0, "negative intron length.");
		cigar[_n_cigar++] = len << 4 | CIGAR_OP_N;
	}


	//last exon
	--iter; p = *iter;

	//len_t = p->dexon->end_t - p->end_t + 1;
	//len_q = p->dexon->end_q - p->end_q + 1;
	len_t = aln->end_t - p->end_t + 1;
	len_q = aln->end_q - p->end_q + 1;
	start_t = p->end_t; start_q = p->end_q;

	ubyte_t* ref_seq = (ubyte_t*)calloc(len_t, 1);
	for (k = start_t, l = 0; k < start_t + len_t && k < l_pac; ++k)
		ref_seq[l++] = pacseq[k>>2] >> ((~k&3)<<1) & 3;

	path_t* path = (path_t*)calloc(len_t+len_q, sizeof(path_t));

	aln_global_core(ref_seq, len_t, (ubyte_t*)(seq+start_q), len_q, &ap, path, &path_len);
	exon_cigar = aln_path2cigar32(path, path_len, &n_exon_cigar);

	if (max_cigar <= _n_cigar + n_exon_cigar) {
		max_cigar = _n_cigar + n_exon_cigar + 1;
		cigar = (uint32_t*)realloc(cigar, sizeof(uint32_t) * max_cigar);
	}

	for (int i = 0; i < n_exon_cigar; ++i)
		cigar[_n_cigar++] = exon_cigar[i];

	free (ref_seq); free (path); free (exon_cigar);

	//cigar[n++] = len << 4 | CIGAR_OP_M;

	*n_cigar = _n_cigar;

	return cigar;
}



jigsaw_cigar_t *jigsaw_spliced_aln2cigar(const jigsaw_spliced_aln_t *aln, ubyte_t *seq,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,
		int band_width, int *n_cigar)
{
	uint32_t *cigar32;
	jigsaw_cigar_t *cigar;
	int i;
	cigar32 = spliced_aln2cigar32(aln, seq, l_pac, pacseq, ntpac, band_width, n_cigar);
	cigar = (jigsaw_cigar_t*)cigar32;
	for (i = 0; i < *n_cigar; ++i)
		cigar[i] = __cigar_create( (cigar32[i]&0xf), (cigar32[i]>>4) );
	return cigar;
}


/*align one query sequence that is spliced
 * n_occ: the number of alignments to report
 */
void jigsaw_aln_one_spliced (bwt_t *const bwt[2], bwa_seq_t *seq, const int *g_log_n,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,
		/*gap_stack_t *stack,*/ int n_occ, const gap_opt_t *opt)
{
	float max_aln_score = 1;
	uint32_t extend_exact = 0, global = 1;

	//int i, n_aln = 0, max_aln = 0;

	seq->c1 = 0;

	if (seq->len < 2 * opt->word_size) return;
	//the size of the query sequence is too small; do nothing

	/*get all words in the query sequence*/
	jigsaw_set_read_words (seq, opt->word_size, opt->word_max_overlap,  opt->mode & BWA_MODE_COMPREAD);

	//fprintf (stderr, "collect words\n");

	/*get all hits of all words in an array*/
	//each hits record the coordinates in query and target sequence, as well as the wid in the query
	list<jigsaw_word_hit_t*> hits;
	jigsaw_collect_word_hits (bwt, seq->words, seq->n_words, g_log_n,
		l_pac, pacseq, ntpac, opt->max_word_occ, opt, &hits);

	int n_hits = hits.size();
	if (n_hits == 0) return;

	jigsaw_sort_word_hits_by_strand_pos_t (&hits);

	//fprintf (stderr, "group words to clusters\n");

	list<jigsaw_spliced_aln_cluster_t *> clusters;
	jigsaw_group_hits_to_spliced_aln_clusters (&hits, opt->max_intron_size + 2 * opt->word_size,
			seq->words, seq->n_words, opt->word_size, l_pac, &clusters);

	//int n_clusters = clusters.size();
	jigsaw_sort_spliced_aln_cluster_by_score (&clusters);

	//get all candidate alignment
	list<jigsaw_spliced_aln_t*> aln;
	list<jigsaw_spliced_aln_cluster_t *>::iterator iter;

	//fprintf (stderr, "search spliced align\n");

	int max_cluster_to_search = 100; //TODO: need to be moved
	int i;
	//TODO: heuristic, if the uniqueness score, or rank of the cluster is too low, stop
	for (i = 0, iter = clusters.begin(); iter != clusters.end() && (*iter)->score <= max_aln_score && i < max_cluster_to_search; ++i, ++iter)
		jigsaw_search_spliced_aln (*iter, seq, opt, extend_exact, global, l_pac, pacseq, ntpac,bwt,g_log_n, &aln);

	//fprintf (stderr, "seqname:%s\tclustsize:%d\talnsize:%d\t", seq->name, int(clusters.size()), int(aln.size()) );
	//recover exonic alignments here, instead of doing it in jigsaw_search_spliced_aln
	for (i = 0, iter = clusters.begin(); iter != clusters.end() && (*iter)->score <= max_aln_score && i < max_cluster_to_search; ++i, ++iter)
	   jigsaw_search_exonic_aln(*iter, seq, opt, &aln) ;
	//fprintf (stderr,"size\t%d\n", int(aln.size()) );
	if (aln.size()> 0) {

		jigsaw_sort_spliced_aln_by_diff (&aln);

		//set the best alignment
		jigsaw_spliced_aln_t *top = aln.front();
		uint32_t best_diff = top->diff;
		//equally best alignments
		int i;
		list<jigsaw_spliced_aln_t*>::iterator aln_iter;
		for (i = 0, aln_iter = aln.begin();
			aln_iter != aln.end () && (*aln_iter)->diff == best_diff
			&& (*aln_iter)->diff <= opt->max_diff
			&& (uint32_t)((*aln_iter)->n_gapo_t + (*aln_iter)->n_gapo_q) <= opt->max_gapo
			&& (uint32_t)((*aln_iter)->n_gapo_t + (*aln_iter)->n_gapo_q + (*aln_iter)->n_gape_t + (*aln_iter)->n_gape_q) <= opt->max_gape
			&& (*aln_iter)->logistic_prob > opt->min_logistic_prob;
			++aln_iter, ++i) ;

		seq->c1 = i;

		//sub optimal alignments
		for (; aln_iter != aln.end ()
			&& (*aln_iter)->diff <= opt->max_diff
			&& (uint32_t)((*aln_iter)->n_gapo_t + (*aln_iter)->n_gapo_q) <= opt->max_gapo
			&& (uint32_t)((*aln_iter)->n_gapo_t + (*aln_iter)->n_gapo_q + (*aln_iter)->n_gape_t + (*aln_iter)->n_gape_q) <= opt->max_gape
			&& (*aln_iter)->logistic_prob > opt->min_logistic_prob;
			++aln_iter, ++i) ;

		seq->c2 = i - seq->c1;
	//	fprintf (stderr,"c1:%u\tc2:%u\n", seq->c1, seq->c2);

		if (seq->c1 > 0) {
			seq->type = seq->c1 == 1 ? BWA_TYPE_UNIQUE : BWA_TYPE_REPEAT;

			//set main alignment
			seq->pos = top->start_t;
			seq->strand = top->strand;
			
			seq->sense_strand = top->sense_strand;

			if(top->junctions != NULL){
			    seq->cigar = jigsaw_spliced_aln2cigar(top, seq->strand ? seq->rseq : seq->seq,
				l_pac, pacseq, ntpac, opt->max_diff, &seq->n_cigar);
			}
			else{
			    seq->cigar = jigsaw_exonic_aln2cigar(top, seq->strand ? seq->rseq : seq->seq,
                                    l_pac, pacseq, ntpac, opt->max_diff, &seq->n_cigar); 
			}

			seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, opt->max_diff, g_log_n);
			//seq->nm = best_diff;
			int nm;
			kstring_t *str = (kstring_t*)calloc(1, sizeof(kstring_t));
			seq->md = bwa_cal_md1(seq->n_cigar, seq->cigar, top->end_q-top->start_q+1, seq->pos, (seq->strand? seq->rseq : seq->seq) + top->start_q,
								l_pac, ntpac? ntpac : pacseq, str, &nm);
			seq->nm = nm;
			//if(nm > (int)opt->max_diff) return;
			//ensure not to output alignments with many mismatches, however, the problem is sometimes aln->diff is over-counted because of overlap between junctions and exons -wuj
			free(str->s); free(str);

			seq->n_mm = top->n_mm;
			seq->n_gapo_t = top->n_gapo_t; seq->n_gapo_q = top->n_gapo_q;
			seq->n_gape_t = top->n_gape_t; seq->n_gape_q = top->n_gape_q;

			//seq->multi = 0; seq->n_multi = 0;
					
		}

		if( i > 1 && 1){ //i = c1 + c2
		    seq->n_multi = i < opt->max_report_multi ? (i - 1) : (opt->max_report_multi -1) ;
		    seq->multi = (bwt_multi1_t*) calloc (seq->n_multi, sizeof(bwt_multi1_t) );
		    //the multi mappers will be recorded here
		    aln_iter = aln.begin();
		    aln_iter++;
		    int nm_tmp;
		    char *md_tmp;
		    kstring_t *str_tmp = (kstring_t*)calloc(1, sizeof(kstring_t));
		    int j=0;
		    for(; j < seq->n_multi; aln_iter++, j++ ) {
			jigsaw_spliced_aln_t *curr_aln = *aln_iter;
			bwt_multi1_t *q = seq->multi + j;
			q->pos = curr_aln->start_t;
			q->strand = curr_aln->strand;
			q->sense_strand = curr_aln->sense_strand;
			int n_cigar_tmp;
			if(curr_aln->junctions != NULL){
			    q->cigar = jigsaw_spliced_aln2cigar(curr_aln, q->strand ? seq->rseq : seq->seq,
				l_pac, pacseq, ntpac, opt->max_diff, &n_cigar_tmp);
			}
			else{
			     q->cigar = jigsaw_exonic_aln2cigar(curr_aln, q->strand ? seq->rseq : seq->seq,
	                        l_pac, pacseq, ntpac, opt->max_diff, &n_cigar_tmp);
			}
			q->n_cigar = n_cigar_tmp;
			md_tmp = bwa_cal_md1(q->n_cigar, q->cigar, curr_aln->end_q - curr_aln->start_q+1, q->pos, (q->strand? seq->rseq : seq->seq) + curr_aln->start_q, l_pac, ntpac? ntpac : pacseq, str_tmp, &nm_tmp);

			q->mm = curr_aln->n_mm;
			q->nm = nm_tmp;// this might be different to curr_aln->diff, due to the global adjustment in the concat step, the same reason curr_aln->diff might be differnt to q->mm+ q->gap_t+ q->gap_q 

			q->gap_t = curr_aln->n_gapo_t + curr_aln->n_gape_t;
			q->gap_q = curr_aln->n_gapo_q + curr_aln->n_gape_q;

			//if ( curr_aln->diff != nm_tmp) fprintf (stderr, "%s\t%u\t%u\t%u\t%d\t%d\n",seq->name, q->mm, q->gap_t, q->gap_q, nm_tmp, curr_aln->diff);
		    }
		    free(str_tmp->s); free(str_tmp);
		    
		}
	}

	jigsaw_destroy_word_hits (&hits);
	jigsaw_destroy_spliced_aln_clusters (&clusters);
	jigsaw_destroy_spliced_aln (&aln);

}



/*Align one query sequence
 * It will first try to align the query not allowing splicing
 * If failed it will then search for spliced alignment
 */
void jigsaw_aln_one (bwt_t *const bwt[2], bwa_seq_t *seq, const int *g_log_n,
		int64_t l_pac, const ubyte_t *pacseq, const ubyte_t *ntpac,
		gap_stack_t *stack, int n_occ, const gap_opt_t *opt)
{
	bwa_seq_t *p = seq;
	//clock_t t;
	//t = clock();
	//do standard search first
	//fprintf (stderr, "standard alignment\n");
	//fprintf (stderr,"%s\n", seq->name);
	gap_opt_t local_opt = *opt;
	if (opt->non_denovo_search && !(opt->splice_site_map) ) { /*local_opt.max_diff = opt.max_diff; */}
	//this only happens when non-denovo without annotation db
	//keep the original max_diff
	else{
	
	    if(local_opt.max_diff >2) local_opt.max_diff = 2;
	    if(opt->splice_site_map ) local_opt.max_diff = 0;
	    //no matter if there is denovo search, we want to recover as many as junctions if the user provide a db
	    //TODO: need to setup an option, if matches can be found with in the diffs, then we don't look for spliced alns, The higher this number, the faster the program is
	}

	jigsaw_cal_sa_reg_gap(bwt, p, stack, &local_opt);

	//if no hits found in standard search, do spliced search
	if (p->n_aln > 0) {
	//if (1) {
		//set main and multi
		bwa_aln2seq_core(p->n_aln, p->aln, p, 1, n_occ);

		//calculate p->pos, p->seQ
		//jigsaw_cal_pac_pos(bwt[0], p, opt->max_diff, opt->fnr, g_log_n);
		jigsaw_cal_pac_pos(bwt[0], p, local_opt.max_diff, local_opt.fnr, g_log_n);

		//now seqs[i].pos and seqs[i].multi[j].pos points to the 5'end on the positive strand of the genome
		//set p->cigar, p->n_cigar, multi[i]->cigar, multi[i]->n_cigar
		//jigsaw_refine_gapped(l_pac, p, pacseq, ntpac, opt->max_diff);
		jigsaw_refine_gapped(l_pac, p, pacseq, ntpac, local_opt.max_diff);
	}
	else {
		//fprintf (stderr, "spliced alignment\n");
		if ( (! opt->non_denovo_search) || opt->splice_site_map) 
		    jigsaw_aln_one_spliced (bwt, seq, g_log_n, l_pac, pacseq, ntpac, n_occ, opt);

		if (seq->c1 == 0) {
			//no match, call original BWA functions
			//set main and multi
			bwa_aln2seq_core(seq->n_aln, seq->aln, seq, 1, n_occ);

			//calculate p->pos, p->seQ, p->mapQ
			//jigsaw_cal_pac_pos(bwt[0], seq, opt->max_diff, opt->fnr, g_log_n);
			jigsaw_cal_pac_pos(bwt[0], seq, local_opt.max_diff, local_opt.fnr, g_log_n);

			//now seqs[i].pos and seqs[i].multi[j].pos points to the 5'end on the positive strand of the genome
			//set p->cigar, p->n_cigar, multi[i]->cigar, multi[i]->n_cigar
			//jigsaw_refine_gapped(l_pac, seq, pacseq, ntpac, opt->max_diff);
			jigsaw_refine_gapped(l_pac, seq, pacseq, ntpac, local_opt.max_diff);
			
		}
	}
//	fprintf(stderr, "%s\t%.2f\t %ld\n", seq->name, (float)(clock() - t),  CLOCKS_PER_SEC);

}


/*Align a batch of query sequence
 * 1. initialize the stack
 * 2. assign each query sequence to a thread
 * 3. do alignment
 */
void jigsaw_aln_batch(int tid, /*thread id, 0 means no multi-threading*/
		bwt_t *const bwt[2], /*bwt[0] is the forward bwt, and bwt[1] is the reverse bwt*/
		int n_seqs, bwa_seq_t *seqs, /*query sequences*/
		const int *g_log_n,
		int64_t l_pac, const ubyte_t *pacseq, /*packed reference sequence*/
		const ubyte_t *ntpac, /*reference genome in color space or null*/
		const gap_opt_t *opt)
{
	int i, /*max_l = 0, */max_len;
	gap_stack_t *stack;
	//const ubyte_t *seq[2];
	gap_opt_t local_opt = *opt;
	//ubyte_t *pacseq, *ntpac = 0;

	int n_occ = 3;

	// initiate priority stack
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

	//core loop
	for (i = 0; i != n_seqs; ++i) {

		//if (opt->verbose && i % 500 == 0) fprintf (stderr, "%d...", i);

		//fprintf (stderr, "%d ...\n", i);
		bwa_seq_t *p = seqs + i;
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
/*
		if (strcmp (p->name, "HWI-ST155_05162111902026#0") == 0)
		{
			fprintf (stderr, "found %s\n", p->name);
		}
*/

#ifdef HAVE_PTHREAD
		if (opt->n_threads > 1) {
			pthread_mutex_lock(&g_seq_lock);
			if (p->tid < 0) { // unassigned
				int j;
				for (j = i; j < n_seqs && j < i + THREAD_BLOCK_SIZE; ++j)
					seqs[j].tid = tid;
			} else if (p->tid != tid) {
				pthread_mutex_unlock(&g_seq_lock);
				continue;
			}
			pthread_mutex_unlock(&g_seq_lock);
		}
#endif

		jigsaw_aln_one (bwt, p, g_log_n, l_pac, pacseq, ntpac, stack, n_occ, &local_opt);
	}
	//free(seed_w[0]); free(seed_w[1]);
	gap_destroy_stack(stack);
	//free(pacseq); free(ntpac);
}


#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt[2];
	int n_seqs;
	bwa_seq_t *seqs;
	const int *g_log_n;
	const ubyte_t *pacseq, *ntpac /*color space*/;
	int64_t l_pac;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	//bwa_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
	jigsaw_aln_batch (d->tid, d->bwt, d->n_seqs, d->seqs, d->g_log_n, d->l_pac, d->pacseq, d->ntpac, d->opt);

	return 0;
}
#endif


bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa)
{
	bwa_seqio_t *ks;
	if (mode & BWA_MODE_BAM) { // open BAM
		int which = 0;
		if (mode & BWA_MODE_BAM_SE) which |= 4;
		if (mode & BWA_MODE_BAM_READ1) which |= 1;
		if (mode & BWA_MODE_BAM_READ2) which |= 2;
		if (which == 0) which = 7; // then read all reads
		ks = bwa_bam_open(fn_fa, which);
	} else ks = bwa_seq_open(fn_fa);
	return ks;
}

/* rgoya: Temporary clone of aln_path2cigar to accomodate for bwa_cigar_t,
__cigar_op and __cigar_len while keeping stdaln stand alone */
jigsaw_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar)
{
	uint32_t *cigar32;
	jigsaw_cigar_t *cigar;
	int i;
	cigar32 = aln_path2cigar32((path_t*) path, path_len, n_cigar);
	cigar = (jigsaw_cigar_t*)cigar32;
	for (i = 0; i < *n_cigar; ++i)
                cigar[i] = __cigar_create( (cigar32[i]&0xf), (cigar32[i]>>4) );
	return cigar;
}


/*The core function to read and align reads, and print output*/
void jigsaw_aln_core(const char *prefix, /*prefix of the genome index*/
		const char *fn_fa, /*fasta or fastq file of query reads*/
		gap_opt_t *opt, /*alignment options*/
		int argc,
		char *argv[])
{
	//variables to calculate SA intervals
	int i, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t,t0;
	//cpu time
	time_t t0_wall;
	//wall time
	bwt_t *bwt[2];
	int g_log_n[256];
	//splice_site_map_t *splice_site_map = 0;

	//variables to convert SA intervals to final alignment
	//int n_occ = 3;
	bntseq_t *bns, *ntbns = 0;
	ubyte_t *pacseq, *ntpac = 0;

	char *bwa_rg_line = 0, *bwa_rg_id = 0;
	t0_wall = time(NULL);
	t0 = clock();
	t = clock();
	// initialization
	{	
		//load packed genomic sequences
		
		if (opt->rg != 0 && bwa_set_rg(opt->rg, &bwa_rg_id, &bwa_rg_line) < 0) {
			fprintf(stderr, "[%s] malformated @RG line: %s\n", __func__, opt->rg);
			exit (1);
		}
		
		bwase_initialize(g_log_n);
	
		//load packed reference sequence
		if (opt->verbose) fprintf (stderr, "loading the packed reference database ...");

		bns = bns_restore(prefix);
		srand48(bns->seed);	

		pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
		rewind(bns->fp_pac);
		fread(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);

		if (!(opt->mode & BWA_MODE_COMPREAD)) { // in color space; initialize ntpac
			ntbns = bwa_open_nt(prefix);
			ntpac = (ubyte_t*)calloc(ntbns->l_pac/4+1, 1);
			rewind(ntbns->fp_pac);
			fread(ntpac, 1, ntbns->l_pac/4 + 1, ntbns->fp_pac);
		}

		if (opt->verbose) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		//load known splice sites
		if (opt->junction_file) {
			if (opt->verbose) fprintf (stderr, "loading the exon junction database ...");
			opt->splice_site_map = build_splice_site_map (bns, opt->junction_file);
			if (opt->verbose) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
			t = clock();
		}
		//load regression model
		if (opt->regression_file) {
		    if (opt->verbose) fprintf (stderr, "loading the regression model ...");
		    load_model_cfg (opt->regression_file);
		    if (opt->verbose) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		    t = clock();

		}
	}

	{ // load BWT
		if (opt->verbose) fprintf (stderr, "loading BWT and suffix array ...");

		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);

		strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
		free(str);

		if (opt->verbose) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();
	}

	//if (opt->verbose) fprintf (stderr, "aligning reads ...\n");

	bwa_print_sam_SQ(bns, bwa_rg_line);
	//print @PG tags
	printf("@PG\tID:OLego\tVN:%s\tCL:",PACKAGE_VERSION );
        for(int i = 0; i < argc; i++ ){
                printf("%s ", argv[i]);
        }
	 printf("\n");
	
	//print information into stderr as well
	if (opt->verbose) {
		fprintf(stderr, "\nOLego: version %s\n", PACKAGE_VERSION);
                fprintf(stderr, "Compiled at %s, %s\n", __TIME__, __DATE__);
	}
	//set ks
	ks = bwa_open_reads(opt->mode, fn_fa);

	// core loop
	///fwrite(opt, sizeof(gap_opt_t), 1, stdout);
	int iter = 0;
	
	//int n_batch = 0x4000;
	//int n_batch = 0x40000;

	while ((seqs = bwa_read_seq(ks, opt->n_batch, &n_seqs, opt->mode & BWA_MODE_COMPREAD, opt->trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();

		if (opt->verbose) fprintf(stderr, "[batch %d]\n", iter);
		iter++;

		if (opt->verbose) fprintf(stderr, "align reads... ");


#ifdef HAVE_PTHREAD
		if (opt->n_threads <= 1) { // no multi-threading at all
			//perform alignment of the batch
			//the result is ready for output
			jigsaw_aln_batch (0, bwt, n_seqs, seqs, g_log_n, bns->l_pac, pacseq, ntpac, opt);
		} else {
			pthread_t *tid;
			pthread_attr_t attr;
			thread_aux_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
			for (j = 0; j < opt->n_threads; ++j) {
				data[j].tid = j; data[j].bwt[0] = bwt[0]; data[j].bwt[1] = bwt[1];
				data[j].n_seqs = n_seqs; data[j].seqs = seqs;
				data[j].g_log_n = g_log_n;
				data[j].l_pac = bns->l_pac;
				data[j].pacseq = pacseq;
				data[j].ntpac = ntpac;
				data[j].opt = opt;
				pthread_create(&tid[j], &attr, worker, data + j);
			}
			for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else
		//bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
		jigsaw_aln_batch (0, bwt, n_seqs, seqs, g_log_n, bns->l_pac, pacseq, ntpac, opt);
#endif

		if (opt->verbose) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		//output
		fprintf(stderr, "print alignments... ");
		for (i = 0; i < n_seqs; ++i)
			bwa_print_sam1(bns, seqs + i, 0, opt->mode, opt->max_top2, bwa_rg_id);

		if (opt->verbose) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		//release memory
		jigsaw_free_read_seq(n_seqs, seqs);
		if (opt->verbose) fprintf(stderr, "%d sequences have been processed.\n", tot_seqs);

	}

	// destroy
	free(pacseq); free(ntpac);
	bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	bwa_seq_close(ks);

	if (ntbns) bns_destroy(ntbns);
	bns_destroy(bns);

	if (opt->splice_site_map) destroy_splice_site_map (opt->splice_site_map);
	if (opt->verbose) fprintf(stderr,"Total cpu time used: \t%.2f sec\n",(float)(clock()-t0) / CLOCKS_PER_SEC);
	if (opt->verbose) fprintf(stderr,"Total wall time used: \t%ld sec\n",(time(NULL)-t0_wall) );
}

