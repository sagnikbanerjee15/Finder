/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef BWT_BNTSEQ_H
#define BWT_BNTSEQ_H

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef uint8_t ubyte_t;
#endif

#define BASE_A 0
#define BASE_C 1
#define BASE_G 2
#define BASE_T 3
#define BASE_N 4

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

extern unsigned char nst_nt4_table[256];

#ifdef __cplusplus
extern "C" {
#endif

	/*inline*/ ubyte_t get_pacseq_base (const ubyte_t *pacseq, int64_t k);
	void bns_dump(const bntseq_t *bns, const char *prefix);
	bntseq_t *bns_restore(const char *prefix);
	bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename);
	void bns_destroy(bntseq_t *bns);
	void bns_fasta2bntseq(gzFile fp_fa, const char *prefix);
	int bns_coor_pac2real(const bntseq_t *bns, int64_t pac_coor, int len, int32_t *real_seq);

#ifdef __cplusplus
}
#endif

#endif
