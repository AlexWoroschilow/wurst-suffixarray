/*
 * 14 June 2005
 * This defines the vector of class membership probability.
 * It lives in its own file since it must be visible to
 * different class or structure functions.
 * It can exist in more than one form. It may be a simple array,
 * or it may be compact, containing only the non-zero
 * probabilities.
 * It may be normalised so that the probabilities sum to 1.0, or
 * so that the sum of squares sums to 1.0.
 * rcsid = $Id: prob_vec.h,v 1.5 2006/02/22 14:00:05 torda Exp $
 */

#ifndef PROB_VEC_H
#define PROB_VEC_H

extern const char PVEC_TRUE_PROB;             /* Normalised to sum = 1.0 */
extern const char PVEC_UNIT_VEC;              /* or to unit vector size */
extern const char PVEC_CRAP;                  /* Don't know yet. */

/*
 * In the compact form cmpct_xxxx,
 * cmpct_n says how many probabilities are used at each site.
 * cmpct_ndx says indices of classes which are used.
 * cmpct_prob is a flat array of corresponding probabilities.
 * 
 * mship[a][b] is the membership of "site" a in class b.
 * The definition of site "a" may be subject to some change in the future.
 * It is, however, dimensioned as mship [n_pvec][n_class].
 *
 * Originally, n_pvec = size - frag_len + 1
 * At some time, we will be able to handle partially described probability
 * vectors, so 
 */
  
struct prob_vec {
    unsigned short int *cmpct_n;  /* Number of stored probabilities per site */
    float *cmpct_prob;                     /* The list of probability values */
    unsigned short int *cmpct_ndx;                   /* Indices within class */
                                                         /* cmpct_n per site */
    float **mship;                  /* The expanded, simple array membership */
    size_t n_pvec;                          /* Number of probability vectors */
    size_t n_class;            /* How many classes are in the classification */
    size_t prot_len;                   /* The length of the original protein */
    size_t frag_len;         /* The length of the fragment. Typically 4 to 9 */
    char norm_type;                     /* True probability or unit vector ? */
    char *compnd;		            /*The compound information string*/
    size_t compnd_len;	             /*The length of the compound information*/
};

#endif /* PROB_VEC_H */
