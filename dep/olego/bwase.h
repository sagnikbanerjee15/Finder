#ifndef BWASE_H
#define BWASE_H

#include "bntseq.h"
#include "bwt.h"
#include "kstring.h"
//#include "bwtaln.h"
#include "jigsawaln.h"

#ifdef __cplusplus
extern "C" {
#endif

	int bwa_set_rg(const char *s, char **bwa_rg_id, char **bwa_rg_line);
	bntseq_t *bwa_open_nt(const char *prefix);

	// Initialize mapping tables in the bwa single-end mapper.
	void bwase_initialize(int *g_log_n);
	// Calculate the approximate position of the sequence from the specified bwt with loaded suffix array.
	void bwa_cal_pac_pos_core(const bwt_t* bwt, bwa_seq_t* seq, const int max_mm, const float fnr, const int *g_log_n);
	// Refine the approximate position of the sequence to an actual placement for the sequence.
	//void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq, bntseq_t *ntbns);
	void jigsaw_refine_gapped(int64_t l_pac, bwa_seq_t *seq, const ubyte_t *pacseq, const ubyte_t *ntpac, int band_width);
	// Backfill certain alignment properties mainly centering around number of matches.

	int bwa_approx_mapQ(const bwa_seq_t *p, int mm, const int *g_log_n);
	int jigsaw_cal_diff (const path_t *path, int path_len,
			const ubyte_t *ref_seq, int ref_len, const ubyte_t *seq, int len,
			uint32_t *n_mm, uint32_t *n_gapo_t, uint32_t *n_gapo_q, uint32_t *n_gape_t, uint32_t *n_gape_q);
	//calculate editting distance of refined alignment


	char *bwa_cal_md1(int n_cigar, jigsaw_cigar_t *cigar, int len, bwtint_t pos, ubyte_t *seq,
					  bwtint_t l_pac, const ubyte_t *pacseq, kstring_t *str, int *_nm);

	void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
	void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);

	// Calculate the end position of a read given a certain sequence.
	int64_t pos_end(const bwa_seq_t *p);

	void jigsaw_cal_pac_pos (const bwt_t *bwt, bwa_seq_t *seq, int max_mm, float fnr, const int *g_log_n);
	int jigsaw_word_hits_total (const jigsaw_word_t *w);
	void jigsaw_collect_word_hits_core (const bwt_t *bwt, jigsaw_word_t *w, list<jigsaw_word_hit_t*> *hits);

	void bwa_print_sam_SQ(const bntseq_t *bns, const char *bwa_rg_line);
	void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2, const char *bwa_rg_id);

	//void bwa_sai2sam_se_core(const char *prefix, const char *fn_sa, const char *fn_fa, int n_occ, const char *bwa_rg);

#ifdef __cplusplus
}
#endif

#endif // BWASE_H
