/*
    This is for calculating the splicesite score for each junction
*/
#include <math.h>
#include "bntseq.h"

// this is for the definition of BASE_X

#include "bwaseqio.h"
// for seq_reverse
/*
#define LOGIT_A -2.388e+00
#define LOGIT_B -4.836e-05
#define LOGIT_C 4.154e-01
#define FOREGROUND_NUM 36262
*/
/* 
#define LOGIT_A 1.127e+00
#define LOGIT_B -4.810e-05
#define LOGIT_C 2.771e-01
#define FOREGROUND_NUM 1
*/
using namespace std;

extern double global_logit_a;
extern double global_logit_b;
extern double global_logit_c;
extern double global_foreground_num;

extern double global_score_matrix[4][60];
// to store matrix of the motif around known splice site, the +-15nt around 5ss and 3ss are connected to present the junction site

extern double global_score_background[4];

double get_seqs_score (int64_t fragment_size, ubyte_t *ref_seq, double (&global_score_matrix)[4][60], double *global_score_background);
// get the score by comparing seq with the matrix of the motif

double get_splice_score(const ubyte_t *pacseq, int64_t ue_end_t, int64_t de_start_t);
// calculate the score using logit regression model

//void load_model_cfg(double *global_logit_a, double *global_logit_b, double *global_logit_c, double *global_foreground_num, double global_score_matrix[4][60], double *global_score_background, char *cfg_file_name);
void load_model_cfg(const char *cfg_file_name);
// load configuration file for the regression model

