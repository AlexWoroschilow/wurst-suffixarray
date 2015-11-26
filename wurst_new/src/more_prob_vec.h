/*
 * 08 Jun 2006
 * rcsid = $Id: empty_prop_vec.h,v 0.1 2006/06/08 16:43:36 ahattesohl Exp $
 */

#ifndef MORE_PROB_VEC_H
#define MORE_PROB_VEC_H

struct aa_strct_clssfcn;
struct pair_set;
struct seq;

struct afloat {
    float a;
};
struct float_array {
    struct afloat *floats;
};

struct prob_vec  *prob_vec_add (struct prob_vec *p_vec1,
                                struct prob_vec *p_vec2,
                                struct pair_set *p_set,
                                const size_t norm_flag);
struct prob_vec  *prob_vec_add2 (struct prob_vec *p_vec1,
                                 struct prob_vec *p_vec2,
                                 struct pair_set *p_set,
                                 const size_t cur_step);
/*
struct prob_vec  *prob_vec_avg (struct prob_vec *p_vec,
                                float **coverage,
                                const size_t norm_flag);
*/
struct prob_vec  *prob_vec_add_weighted (struct prob_vec *p_vec1,
                                         struct prob_vec *p_vec2,
                                         struct pair_set *p_set,
                                         const size_t norm_flag);
void              prob_vec_2_true_prob (struct prob_vec *p_v);
void              prob_vec_2_unit_vec (struct prob_vec *p_v);
float             pair_set_get_seq_id (struct pair_set *p_set,
                                       struct seq *seq1,
                                       struct seq *seq2);

#endif /*  MORE_PROB_VEC_H */
