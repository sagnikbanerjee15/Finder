#ifndef BWASEQIO_H
#define BWASEQIO_H

#include "bntseq.h"
#include "bwt.h"
//#include "bwtaln.h"
#include "jigsawaln.h"

#ifdef __cplusplus
extern "C" {
#endif

	bwa_seqio_t *bwa_bam_open(const char *fn, int which);
	bwa_seqio_t *bwa_seq_open(const char *fn);
	void bwa_seq_close(bwa_seqio_t *bs);
	bwa_seq_t *bwa_read_seq(bwa_seqio_t *bs, int n_needed, int *n, int is_comp, int trim_qual);
	//static bwa_seq_t *bwa_read_bam(bwa_seqio_t *bs, int n_needed, int *n, int is_comp, int trim_qual);

	void seq_reverse(int len, ubyte_t *seq, int is_comp);
	//static int bwa_trim_read(int trim_qual, bwa_seq_t *p);

	void jigsaw_set_read_words (bwa_seq_t *p, int word_size, int word_max_overlap,  int is_comp);
	//static void jigsaw_free_read_seq_one (bwa_seq_t *p);
	void jigsaw_free_read_seq(int n_seqs, bwa_seq_t *seqs);


#ifdef __cplusplus
}
#endif

#endif // BWASEQIO_H



