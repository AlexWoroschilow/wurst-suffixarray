/*
 * 25 Sep 2001
 * Some operations on pair_set's (alignments)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "e_malloc.h"
#include "mprintf.h"
#include "pair_set_i.h"
#include "pair_set.h"
#include "read_seq_i.h"
#include "scratch.h"
#include "seq.h"
#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: pair_set.c,v 1.20 2007/06/01 12:22:18 mmosisch Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */

/* ---------------- pair_set_string    --------------------------
 * This could move to another file. It needs to know about
 * the innards of the sequence structure. Really, this file
 * should only contain stuff for pair sets, regardless of their
 * type.
 */
char *
pair_set_string (struct pair_set *pair_set, struct seq *s1, struct seq *s2)
{
    int **indices = pair_set -> indices;
    size_t i;
    size_t j;
    const char GAP = '-';

    seq_thomas2std (s1);
    seq_thomas2std (s2);

    scr_reset();
    for (j = 0; j < pair_set->m; j++) {
        for (i = 0; i < pair_set->n; i++) {
            char c1 = (indices[i][j] == GAP_INDEX ? GAP : 'X'); /*s1->seq[indices[i][j]]);*/
            scr_printf ("%c", c1);
        }
        scr_printf ("\n");
    }

    return (scr_printf ("%c", '\n'));
}

/* ---------------- multal_string    --------------------------
 * This could move to another file. It needs to know about
 * the innards of the sequence structure. Really, this file
 * should only contain stuff for pair sets, regardless of their
 * type.
 */
char *
multal_string (struct pair_set *pair_set)
{
    int **indices = pair_set -> indices;
    size_t i;
    size_t j;
    const char GAP = '-';

    scr_reset();
    for (j = 0; j < pair_set->m; j++) {
		for (i = 0; i<indices[0][j]; i++){
			scr_printf ("%c", GAP);
		}
        for (i = 0; i < pair_set->n; i++) {
            char c1 = (indices[i][j] == GAP_INDEX ? GAP : 'X'); /*s1->seq[indices[i][j]]);*/
            scr_printf ("%c", c1);
        }
        scr_printf ("\n");
    }

    return (scr_printf ("%c", '\n'));
}

/* ---------------- pair_set_score    --------------------------
 */
int
pair_set_score (struct pair_set *s, float *score, float *scr_smpl)
{
    if (s == NULL)
        return EXIT_FAILURE;
    *score    = s->score;
    *scr_smpl = s->smpl_score;
    return EXIT_SUCCESS;
}

/* ---------------- pair_set_extend   --------------------------
 * Take a short alignment and extend it. Well documented in
 * in file wurst.pod.
 */
int
pair_set_extend (struct pair_set *s, const size_t n0,
                 const size_t n1, const int long_or_short)
{
    int **p, *plast, **newpairs, **np;
    const char *this_sub = "pair_set_extend";
    size_t left, right, long_left, long_right, newsize, idx;
    int long_left_on_a  = -999,   /* The initialisation is only to provoke */
        long_right_on_a = -999;   /* a crash if code below is not correct */

    long_left = long_right = 0; /* Not necessary, but stops compiler warning */

    if ((long_or_short != EXT_SHORT) && (long_or_short != EXT_LONG)) {
        err_printf (this_sub, "Must be fed either $EXT_LONG or $EXT_SHORT\n");
        return EXIT_FAILURE;
    }
    if (s->n == 0) {
        if (long_or_short == EXT_LONG) {
            size_t i;
            newsize = n0 + n1;
            newpairs = i_matrix(newsize, 2);
            for (idx = 0, i = 0; i < n0; i++, idx++) {
                np[idx][0] = i;
                np[idx][1] = GAP_INDEX;
            }
            for (i = 0; i < n1; i++, np+=2) {
                np[idx][0] = GAP_INDEX;
                np[idx][1] = i;
            }
        } else {
            newsize = 0;
            newpairs = NULL;
        }
        goto go_back;
    }
    p = s->indices;
    plast = s->indices[s->n - 1];      /* last of pairlist */

    if (p[0][0] < 0 || p[0][1] < 0 || p[s->n - 1][0] < 0 || p[s->n - 1][1] < 0) {
        err_printf (this_sub, "This pair set has already been extended\n");
        return EXIT_FAILURE;
    }

    {   /* for short or long extension, we need left and right values. */
        size_t ia, ib;
        left = (p[0][0] < p[0][1]) ? p[0][0] : p[0][1];
        ia = n0 - p[s->n - 1][0];
        ib = n1 - p[s->n - 1][1];
        right = ((ia < ib) ? ia : ib) - 1;
    }
    newsize = s->n + left + right;
    if (long_or_short == EXT_LONG) {
        size_t a_dangle, b_dangle;
        if (p[0][0] > p[0][1]) {
            long_left = p[0][0];
            long_left_on_a = 1;
        } else {
            long_left = p[0][1];
            long_left_on_a = 0;
        }
        long_left -= left;
        a_dangle = n0 - p[s->n - 1][0];
        b_dangle = n1 - p[s->n - 1][1];
        if (a_dangle > b_dangle) {
            long_right = n0 - p[s->n - 1][0];
            long_right_on_a = 1;
        } else {
            long_right = n1 - p[s->n - 1][1];
            long_right_on_a = 0;
        }
        long_right = long_right - right - 1;
        newsize += long_left + long_right;
    }

    newpairs = i_matrix(newsize, 2);
    np = newpairs;
    if ((long_or_short == EXT_LONG) && (long_left > 0)) { /* long left ext'n */
        size_t i;
        if (long_left_on_a) {
            for (i = 0, idx=0; i < long_left; i++, idx++) {
                np[idx][0] = i;
                np[idx][1] = GAP_INDEX;
            }
        } else {
            for (i = 0, idx=0; i < long_left; i++, idx++) {
                np[idx][0] = GAP_INDEX;
                np[idx][1] = i;
            }
        }
    }

    {                                                /* short left extension */
        int ia = p[0][0] - left;
        int ib = p[0][1] - left;
        for (idx = 0 ; ia < p[0][0]; ia++, ib++, idx++) {
            np[idx][0] = ia;
            np[idx][1] = ib;
        }
    }
	kill_i_matrix(s->indices);
    /*memcpy (np, s->indices, s->n * sizeof (np[0][0]) * 2);*/
	s->indices = copy_i_matrix(np, newsize, 2);
    /*np += s->n;*/
    idx += s->n;
    if (right > 0) {                                /* short right extension */
        size_t i;
        int ia = p[s->n][0] + 1;
        int ib = p[s->n][1] + 1;
        for ( i = 0; i < right ; i++, idx++, ia++, ib++) {
            np[idx][0] = ia;
            np[idx][1] = ib;
        }
    }

    if (long_or_short == EXT_LONG && long_right > 0) {/* long right extension*/
        if (long_right_on_a) {
            size_t i = ((np - 2)[0][0]) + 1;
            for ( ; i < n0; idx++, i++) {
                np[idx][0] = i;
                np[idx][1] = GAP_INDEX;
            }
        } else {
            size_t i = np[idx-2][1] + 1;
            for ( ; i < n1; idx++, i++) {
                np[idx][0] = GAP_INDEX;
                np[idx][1] = i;
            }
        }
    }
    kill_i_matrix (s->indices);
 go_back:
    s->indices = newpairs;
    s->n = newsize;
    return EXIT_SUCCESS;
}

/* ---------------- pair_set_coverage --------------------------
 * Given a pair_set structure and a sequence, return a string
 * which tells us which sites in the sequence and structure are filled.
 * We allocate the strings with E_MALLOC and return a pointer to
 * each char * via pcover1 and pcover2. The caller knows that is
 * has to free() the memory later.
 * Maybe we have done a sequence to structure alignment, maybe something
 * else. Whatever the case, we might be asked about coverage from
 * an alignment.
 */
int
pair_set_coverage (struct pair_set *p_s, size_t s1, size_t s2,
                   char **pcover1, char **pcover2)
{
    char *s_a, *s_b;
    int **p;
    size_t i;
    size_t size_s_a = sizeof (s_a[0]) * (s1 + 1);
    size_t size_s_b = sizeof (s_b[0]) * (s2 + 1);
    const char *this_sub = "pair_set_coverage";
    const char *too_small =
        "The sizes, %d and %d are too small for the pair_set.\n";

    s_a = E_MALLOC (size_s_a);
    memset (s_a, '0', size_s_a); /* Default (not covered) is '0' */
    s_a[s1] = '\0';              /* Will end up as string in interp, so */
                                 /* needs NULL terminator */
    s_b = E_MALLOC (size_s_b);
    memset (s_b, '0', size_s_b);
    s_b[s2] = '\0';

    p = p_s->indices;
    /*plast = p + p_s->n*2;   Number of pairs in list */
    for ( i=0 ; i < p_s->n; i++ ) {
        if (p[i][0] == GAP_INDEX)
            continue;
        if (p[i][1] == GAP_INDEX)
            continue;
        if (p[i][0] > (int) s1)
            goto error;
        if (p[i][1] > (int) s2)
            goto error;
        s_a[p[i][0]] = '1';
        s_b[p[i][1]] = '1';
    }
    *pcover1 = s_a;
    *pcover2 = s_b;
    return EXIT_SUCCESS;
 error:
    s_a = E_REALLOC (s_a, 1);
    s_b = E_REALLOC (s_a, 1);
    s_a[0] = s_b[0] = '\0';
    err_printf (this_sub, too_small, s1, s2);
    err_printf (this_sub, "Seriously broken\n");
    return EXIT_FAILURE;
}

/* ---------------- pair_set_gap      --------------------------
 * Calculate gaps in a sequence, for use by model rescoring
 * code. We return two costs:
 *  gap opening
 *  gap widening.
 * To do this, we need two scale values.
 * We return each of the costs separately. If we were not
 * interested in optimising parameters, this routine should
 * combine the penalties and return one number.
 */
int
pair_set_gap (struct pair_set *p_s, float *open_cost, float *widen_cost,
              const float open_scale, const float widen_scale)
{
    int **p;
    unsigned open, widen;
    size_t idx;
    enum {ALIGNED, INGAP} state;
    const char *this_sub = "pair_set_gap";
    open = 0;
    widen = 0;

    if (p_s == NULL) {
        err_printf (this_sub, "pair_set broken\n");
        return EXIT_FAILURE;
    }
    p = p_s->indices;
    /*plast = p + p_s->n*2;*/
    state = ALIGNED;
    for ( idx = 0; idx < p_s->n; idx++) {
        switch (state) {
        case ALIGNED:
            if (p[idx][1] == GAP_INDEX) {
                state = INGAP;
                open++;
            }
            break;
        case INGAP:
            if (p[idx][1] == GAP_INDEX)
                widen++;
            else
                state = ALIGNED;
        }
    }

    *open_cost  = open * open_scale;
    *widen_cost = widen * widen_scale;

    return EXIT_SUCCESS;
}

/* ---------------- pair_set_destroy  --------------------------
 */
void
pair_set_destroy (struct pair_set *p_s)
{
    if (p_s == NULL)
        return;
    if ( p_s->indices ){
        kill_i_matrix( p_s->indices );
    }
    free (p_s);
}


/* ----------------- pair_set_circularpermutated ------------
 * 0 false
 * 1 true
 */
 int
 pair_set_circularpermutated( struct pair_set *p_s, size_t prot1_len, size_t prot2_len ){

    int bpermutated = 0; /*boolean*/
	int first_boundary  = (int) prot1_len - 1;
	int second_boundary = (int) prot2_len - 1;
    int prot1_start_index, prot1_stop_index, prot2_start_index, prot2_stop_index ;

    if ( p_s->n > 0 )
    {
    	prot1_start_index = p_s->indices[ 0 ][ 0 ];
		prot1_stop_index  = p_s->indices[ p_s->n - 1 ][ 0 ];
		prot2_start_index = p_s->indices[ 0 ][ 1 ];
		prot2_stop_index  = p_s->indices[ p_s->n - 1 ][ 1 ];

		/* proof, if the alignment goes over the trivial matrix-size   */
		if ( ((prot1_start_index <= first_boundary ) && (prot1_stop_index > first_boundary))
		 ||  ((prot2_start_index <= second_boundary) && (prot2_stop_index > second_boundary)) ){
		 	bpermutated = 1; /* true */
		 }

		 /* proof, if any region is aligned doubled,
		  * if the value of: (startpos - stoppos) * (-1) > startpos, then
		  * something will be doubled aligned and is not permitted !!
		  */

	     if (
	          ( (prot1_start_index - prot1_stop_index) * (-1) < prot1_start_index )
	      ||  ( (prot2_start_index - prot2_stop_index) * (-1) < prot2_start_index ) ){
	      	 bpermutated = 0; /* false */
	      }
 	}

	return bpermutated;
 }

/* ---------------- pair_set_get_alignment_indices   --------------------------
 * To give an idea of the size of elements in a score matrix.
 * Do not look at first or last rows or columns.
 */
void
pair_set_get_alignment_indices (struct pair_set *p_s, int sequencenumber, int *start, int *stop)
{
    static const char *this_sub = "pair_set_get_alignment";
    if ( (sequencenumber <= (int) p_s->m) && (p_s->n > 0 ) ){
        *start = p_s->indices[  0  ][ sequencenumber ];
        *stop  = p_s->indices[ (int) (p_s->n - 1)  ][ sequencenumber ];
    } else {
        *start = *stop = 0;
        err_printf (this_sub, "sequencenumber too high \n");
    }
}
