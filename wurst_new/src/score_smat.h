/*
 * 15 Sep 2001
 * rcsid = $Id: score_smat.h,v 1.6 2004/02/26 11:40:20 torda Exp $
 */
#ifndef SCORE_SMAT_H
#define SCORE_SMAT_H

struct score_mat;
struct seq;
struct sub_mat;
struct seqprof;
int
score_smat (struct score_mat *score_mat, struct seq *s1, struct seq *s2,
             const struct sub_mat *smat);
int
score_sprof (struct score_mat *score_mat, struct seqprof *sp,
             struct seq *seq, const struct sub_mat *smat);
int
score_prof_prof (struct score_mat *score_mat, struct seqprof *sp1,
                 struct seqprof *sp2, const struct sub_mat *smat);
#endif /*  SCORE_SMAT_H */
