/*
 * 08 June 2006
 * Functions for handling class membership probability vectors.
 * This could be implemented so the vectors are hidden from the
 * perl/interface level, but it is likely we will want to be able
 * to manipulate them when doing things like sequence
 * optimisation. Thus, we make them fully fledged objects visible
 * to the interpreter.
 */

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>

#include "gsldir/gsl_sf_erf.h"
#include "gsldir/gsl_linalg.h"
#include "gsldir/gsl_blas.h"
#include "amino_a.h"
#include "bad_angle.h"
#include "coord.h"
#include "coord_i.h"
#include "classifyStructure.h"
#include "e_malloc.h"
#include "fio.h"
#include "matrix.h"
#include "mprintf.h"
#include "prob_vec.h"
#include "prob_vec_i.h"
#include "read_ac_strct.h"
#include "read_seq_i.h"
#include "seq.h"
#include "seqprof.h"
#include "str.h"
#include "more_prob_vec.h"
#include "pair_set_i.h"
#include "pair_set.h"


#ifdef temp_deleted_struct_def
struct aa_strct_clssfcn {
    enum {
        SINGLE_NORMAL,   /* gaussian normal distribution            */
        MULTI_NORMAL,    /* correlated gaussian normal distribution */
        UNKNOWN          /* unknown distribution                    */
    } **classmodel;      /* which kinds of models are used          */

    struct oneclass *classes;
    size_t n_class;
    size_t n_case;       /* How much data we had, big number like 750000 */
    size_t n_angle;
    size_t n_aa;         /* number of amino acid descriptors */
    float abs_error;     /* absolut error in (input) measurement */
};
#endif /*  temp_deleted_struct_def */

#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: more_prob_vec.c,v 0.1 2006/06/09 14:11:57 ahattesohl Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */


#ifdef want_empty_prob_vec
struct prob_vec *
empty_prob_vec(size_t size, const struct aa_strct_clssfcn *cmodel)
{
    float *mship_aa, *mship_strct;
    struct prob_vec *pvec;
    const size_t dim = cmodel->n_angle * 2;
    size_t tomall;
    const size_t frag_len   = (size_t) dim / 2;
    const size_t c_size = size;
    const size_t n_pvec = c_size - frag_len + 1;

    tomall = cmodel->n_class * sizeof(mship_aa);
    mship_aa    = E_MALLOC(tomall);
    mship_strct = E_MALLOC(tomall);

    pvec = new_pvec(frag_len, c_size, n_pvec, cmodel->n_class);
    memset(pvec->mship[0], 0, n_pvec * cmodel->n_class * sizeof(mship_aa));
    pvec->norm_type = PVEC_TRUE_PROB;

    free (mship_aa);
    free (mship_strct);

    return pvec;
}
#endif /*  want_empty_prob_vec */


 /* deprecated!!!! */

struct prob_vec *
prob_vec_add(struct prob_vec *p_vec1,
             struct prob_vec *p_vec2,
             struct pair_set *p_set,
             const size_t norm_flag)
{
    struct prob_vec *pvec;
    struct apair *pairs = p_set->pairs;
    size_t frag_len = p_vec1->frag_len;
    size_t prot_len = p_vec1->prot_len;
    size_t n_pvec = p_vec1->n_pvec;
    size_t n_class = p_vec1->n_class;
    size_t align_len = p_set->n;
    float weight1, weight2;
    size_t i, j, k;
    const char *this_sub = "prob_vec_add";

    if (frag_len != p_vec2->frag_len) {
        err_printf (this_sub, "can not add probability vectors with ");
        err_printf (this_sub, "different fragment lengths!\n");
    }
    if (n_class != p_vec2->n_class) {
        err_printf (this_sub, "can not add probability vectors with ");
        err_printf (this_sub, "different number of classes!\n");
    }

    prob_vec_unit_vec (p_vec1);
    prob_vec_unit_vec (p_vec2);
    pvec = new_pvec (frag_len, prot_len, n_pvec, n_class);
    memset (pvec->mship[0], 0, n_pvec * n_class * sizeof (float));
    j = 0;

    if (prob_vec_expand (p_vec1) == EXIT_FAILURE) {
        err_printf (this_sub, "fail on vec 1\n");
        return NULL;
    }
    if (prob_vec_expand (p_vec2) == EXIT_FAILURE) {
        err_printf (this_sub, "fail on vec 2\n");
        return NULL;
    }
    for (i = 0; i < align_len; i++) {
        if (pairs[i].a != GAP_INDEX) {
            if ((pairs[i].b == GAP_INDEX) ||
               ((unsigned) pairs[i].b >= p_vec2->n_pvec)) {
                if (!((unsigned) pairs[i].a >= p_vec1->n_pvec)) {
                    for (k = 0; k < pvec->n_class; k++) {
                        pvec->mship[j][k] =
                        p_vec1->mship[pairs[i].a][k];
                    }
                }
            } else {
                if (!((unsigned) pairs[i].a >= p_vec1->n_pvec)) {
                    weight1 = weight2 = 0.5;
                    for (k = 0; k < pvec->n_class; k++) {
                        pvec->mship[j][k] =
                        (p_vec1->mship[pairs[i].a][k] * weight1) +
                        (p_vec2->mship[pairs[i].b][k] * weight2);
                    }
                }
            }
            j++;
        }
    }

    pvec->cmpct_n = NULL;
    pvec->cmpct_prob = NULL;
    pvec->cmpct_ndx = NULL;

    if (norm_flag) {
        pvec->norm_type = PVEC_TRUE_PROB;
        prob_vec_unit_vec (pvec);
    } else {
        pvec->norm_type = PVEC_CRAP;
    }

    return (pvec);
}


struct prob_vec *
prob_vec_add2(struct prob_vec *p_vec1,
             struct prob_vec *p_vec2,
             struct pair_set *p_set,
             const size_t cur_step)
{
    struct prob_vec *pvec;
    struct apair *pairs = p_set->pairs;
    size_t frag_len = p_vec1->frag_len;
    size_t prot_len = p_vec1->prot_len;
    size_t n_pvec = p_vec1->n_pvec;
    size_t n_class = p_vec1->n_class;
    size_t align_len = p_set->n;
    size_t i, j, k;
    float weight = cur_step / (cur_step + 1);
    const char *this_sub = "prob_vec_add2";

    if (frag_len != p_vec2->frag_len) {
        err_printf (this_sub, "can not add probability vectors with ");
        err_printf (this_sub, "different fragment lengths!\n");
    }
    if (n_class != p_vec2->n_class) {
        err_printf (this_sub, "can not add probability vectors with ");
        err_printf (this_sub, "different number of classes!\n");
    }

    prob_vec_unit_vec (p_vec1);
    prob_vec_unit_vec (p_vec2);

    pvec = new_pvec (frag_len, prot_len, n_pvec, n_class);
    memset (pvec->mship[0], 0, n_pvec * n_class * sizeof (float));
    j = 0;

    if (prob_vec_expand (p_vec1) == EXIT_FAILURE) {
        err_printf (this_sub, "fail on vec 1\n");
        return NULL;
    }
    if (prob_vec_expand (p_vec2) == EXIT_FAILURE) {
        err_printf (this_sub, "fail on vec 2\n");
        return NULL;
    }
    for (i = 0; i < align_len; i++) {
        if (pairs[i].a != GAP_INDEX) {
            if ((pairs[i].b == GAP_INDEX) ||
               ((unsigned) pairs[i].b >= p_vec2->n_pvec)) {
                if (!((unsigned) pairs[i].a >= p_vec1->n_pvec)) {
                    for (k = 0; k < pvec->n_class; k++) {
                        pvec->mship[j][k] = p_vec1->mship[pairs[i].a][k];
                    }
                }
            } else {
                if (!((unsigned) pairs[i].a >= p_vec1->n_pvec)) {
                    for (k = 0; k < pvec->n_class; k++) {
                        pvec->mship[j][k] =
                        (p_vec1->mship[pairs[i].a][k]) * weight +
                        (p_vec2->mship[pairs[i].b][k]) * (1 - weight);
                    }
                }
            }
            j++;
        }
    }

    pvec->cmpct_n = NULL;
    pvec->cmpct_prob = NULL;
    pvec->cmpct_ndx = NULL;

    pvec->norm_type = PVEC_TRUE_PROB;
    prob_vec_unit_vec (pvec);

    return (pvec);
}

/*
struct prob_vec *
prob_vec_avg(struct prob_vec *p_vec,
             float **coverage,
             const size_t norm_flag)
{
    struct prob_vec *pvec;
    size_t frag_len = p_vec->frag_len;
    size_t prot_len = p_vec->prot_len;
    size_t n_pvec = p_vec->n_pvec;
    size_t n_class = p_vec->n_class;
    size_t i, j;
    const char *this_sub = "prob_vec_avg";

    if (norm_flag) {
        prob_vec_unit_vec (p_vec);
    }
    pvec = new_pvec (frag_len, prot_len, n_pvec, n_class);
    memset (pvec->mship[0], 0, n_pvec * n_class * sizeof (float));

    if (prob_vec_expand (p_vec) == EXIT_FAILURE) {
        err_printf (this_sub, "fail on vec\n");
        return NULL;
    }
    for (i = 0; i < n_pvec; i++) {
        for (j = 0; j < pvec->n_class; j++) {
            pvec->mship[i][j] = p_vec->mship[i][j] * *coverage[i];
        }
    }

    pvec->cmpct_n = NULL;
    pvec->cmpct_prob = NULL;
    pvec->cmpct_ndx = NULL;

    pvec->norm_type = PVEC_TRUE_PROB;
    prob_vec_unit_vec (pvec);

    return (pvec);
}
*/

struct prob_vec *
prob_vec_add_weighted(struct prob_vec *p_vec1,
                      struct prob_vec *p_vec2,
                      struct pair_set *p_set,
                      const size_t norm_flag)
{
    struct prob_vec *pvec;
    struct apair *pairs = p_set->pairs;
    size_t frag_len = p_vec1->frag_len;
    size_t prot_len = p_vec1->prot_len;
    size_t n_pvec = p_vec1->n_pvec;
    size_t n_class = p_vec1->n_class;
    size_t align_len = p_set->n;
    float weight1, weight2;
    float weight_step = 1 / (2 * (float) frag_len);
    size_t i, j, k;
    const char *this_sub = "prob_vec_add";

    if (frag_len != p_vec2->frag_len) {
        err_printf (this_sub, "can not add probability vectors with ");
        err_printf (this_sub, "different fragment lengths!\n");
    }
    if (n_class != p_vec2->n_class) {
        err_printf (this_sub, "can not add probability vectors with ");
        err_printf (this_sub, "different number of classes!\n");
    }

    prob_vec_unit_vec (p_vec1);
    prob_vec_unit_vec (p_vec2);
    pvec = new_pvec (frag_len, prot_len, n_pvec, n_class);
    memset (pvec->mship[0], 0, n_pvec * n_class * sizeof (float));
    j = 0;

    if (prob_vec_expand (p_vec1) == EXIT_FAILURE) {
        err_printf (this_sub, "fail on vec 1\n");
        return NULL;
    }
    if (prob_vec_expand (p_vec2) == EXIT_FAILURE) {
        err_printf (this_sub, "fail on vec 2\n");
        return NULL;
    }
    for (i = 0; i < align_len; i++) {
        if (pairs[i].a != GAP_INDEX) {
            if ((pairs[i].b == GAP_INDEX) ||
               ((unsigned) pairs[i].b >= p_vec2->n_pvec)) {
                if (!((unsigned) pairs[i].a >= p_vec1->n_pvec)) {
                    for (k = 0; k < pvec->n_class; k++) {
                        pvec->mship[j][k] =
                        p_vec1->mship[pairs[i].a][k];
                    }
                }
            } else {
                if (!((unsigned) pairs[i].a >= p_vec1->n_pvec)) {
                    weight1 = weight2 = 0.5;

                    for (k = 1; k < frag_len; k++) {
                        if (!(i + k > n_pvec)) {
                            if (pairs[i + k].a == GAP_INDEX) {
                                weight1 -= weight_step;
                                weight2 += weight_step;
                            }
                            if (pairs[i + k].b == GAP_INDEX) {
                                weight1 += weight_step;
                                weight2 -= weight_step;
                            }
                        }
                    }
                    
                    for (k = 0; k < pvec->n_class; k++) {
                        pvec->mship[j][k] =
                        (p_vec1->mship[pairs[i].a][k] * weight1) +
                        (p_vec2->mship[pairs[i].b][k] * weight2);
                    }
                }
            }
            j++;
        }
    }

    pvec->cmpct_n = NULL;
    pvec->cmpct_prob = NULL;
    pvec->cmpct_ndx = NULL;

    if (norm_flag) {
        pvec->norm_type = PVEC_TRUE_PROB;
        prob_vec_unit_vec (pvec);
    } else {
        pvec->norm_type = PVEC_CRAP;
    }

    return (pvec);
}


/* ---------------- prob_vec_true_prob ------------------------
 * Look at each of the probability vectors and normalise the
 * sum to 1.0. This is really two functions, depending on whether
 * we are in compressed or expanded format.
 * This function has been tested rather thoroughly, but is
 * commented out since the current applications do not need it.
 */
void
prob_vec_2_true_prob (struct prob_vec *p_v)
{
    if (p_v -> norm_type == PVEC_TRUE_PROB)
        return;
    if ( p_v->mship) {                           /* We are in expanded form */
        float **mship = p_v->mship;
        float **mlast = mship + p_v->n_pvec;
        for ( ; mship < mlast; mship++) {
            double total = 0.0;
            float *p;
            float *f = *mship;                     /* Point to start of row */
            float *flast = f + p_v->n_class;
            for ( ; f < flast; f++)
                total += *f;
            if (total == 0.0)         /* If broken, then don't touch vector */
                total = 1.;
            for ( p = *mship; p < flast; p++)
                *p /= (float) total;
        }
    } else {                                                /* Compact form */
        float *p_prob, *prob2, *pr_last;
        unsigned short *pi, *p_ndx, *pi_last;
        pi      = p_v->cmpct_n;
        p_ndx   = p_v->cmpct_ndx;
        pi      = p_v->cmpct_n;
        pi_last = pi + p_v->n_pvec;
        p_prob = p_v->cmpct_prob;
        prob2  = p_prob;
        for ( ; pi < pi_last; pi++) {              /* Over each prob vector */
            double  total = 0.0;
            pr_last = p_prob + *pi;
            for ( ; p_prob < pr_last; p_prob++)  /* Each element in one vec */
                total += *p_prob;
            if (total == 0.0)         /* If broken, then don't touch vector */
                total = 1.;
            for ( ; prob2 < pr_last; prob2++)
                *prob2 = (float) (*prob2 / total);
        }
    }
    p_v->norm_type = PVEC_TRUE_PROB;
}


/* ---------------- prob_vec_unit_vec -------------------------
 * If our vectors are normalised so they sum to 1, make the sum
 * of squares by 1.0. That is, normalise the vector to unit
 * length.
 */
void
prob_vec_2_unit_vec (struct prob_vec *p_v)
{
    if (p_v -> norm_type == PVEC_UNIT_VEC)
        return;
    if ( p_v->mship) {                           /* We are in expanded form */
        float **mship = p_v->mship;
        float **mlast = mship + p_v->n_pvec;
        for ( ; mship < mlast; mship++) {
            double total = 0.0;
            float *p;
            float *f = *mship;                     /* Point to start of row */
            float *flast = f + p_v->n_class;
            for ( ; f < flast; f++)
                total += *f * *f;
            if (total == 0.0)         /* If broken, then don't touch vector */
                total = 1;
            total = sqrt (total);
            for ( p = *mship; p < flast; p++)
                *p /= (float) total;
        }
    } else {                                                /* Compact form */
        float *p_prob, *prob2, *pr_last;
        unsigned short *pi, *p_ndx, *pi_last;
        pi      = p_v->cmpct_n;
        p_ndx   = p_v->cmpct_ndx;
        pi      = p_v->cmpct_n;
        pi_last = pi + p_v->n_pvec;
        p_prob = p_v->cmpct_prob;
        prob2  = p_prob;
        for ( ; pi < pi_last; pi++) {              /* Over each prob vector */
            double  total = 0.0;
            pr_last = p_prob + *pi;
            for ( ; p_prob < pr_last; p_prob++)  /* Each element in one vec */
                total += (*p_prob * *p_prob);
            if (total == 0.0)         /* If broken, then don't touch vector */
                total = 1;
            total = sqrt (total);
            for ( ; prob2 < pr_last; prob2++)
                *prob2 = (float) (*prob2 / total);
        }
    }
    p_v->norm_type = PVEC_UNIT_VEC;
}


struct triplet {
    unsigned length;
    unsigned aligned;
    unsigned ident;
};

static struct triplet
get_seq_id ( struct pair_set *pair_set, struct seq *s1, struct seq *s2)
{
    struct triplet result;
    unsigned length, aligned, ident;
    struct apair *p, *plast;
    p = pair_set->pairs;

    length = aligned = ident = 0;
    plast = p + pair_set->n;
    for ( ; p < plast; p++) {
        int a = p->a, b = p->b;
        length++;
        if ((a == GAP_INDEX) || (b == GAP_INDEX))
            continue;
        aligned++;
        if (tolower(s1->seq[a]) == tolower(s2->seq[b]))
            ident++;
    }
    result.length = length;
    result.ident = ident;
    result.aligned = aligned;
    if (length != pair_set->n)
        err_printf ("get_seq_id", "Silly bloody bug %s\n", __FILE__);
    return (result);
}

float
pair_set_get_seq_id (struct pair_set *p_set,
                     struct seq *seq1,
                     struct seq *seq2)
{
    struct triplet seq_info;
    float result;

    seq_info = get_seq_id (p_set, seq1, seq2);
    
    mprintf("length: %u\n", seq_info.length);
    mprintf("ident: %u\n", seq_info.ident);
    mprintf("aligned: %u\n", seq_info.aligned);
    
    result = ((float) seq_info.ident / (float) seq_info.length) * 100;
    mprintf("result: %f\n", result);
    
    return (result);
}
