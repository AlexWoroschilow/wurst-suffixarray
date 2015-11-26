/* 23-1-2004
   functions operating with
   a scor_set
 */

#ifndef SCOR_SET_I_H
#define SCOR_SET_I_H

struct scor_set;
struct pair_set;
struct score_mat;
struct scor_set *scor_set_simpl (struct pair_set *pair_set, const struct score_mat *smat);
struct scor_set *scor_set_fromvec ( size_t n, double *v );
int              scor_set_scale(struct scor_set *ss, float scale);
void             scor_set_destroy ( struct scor_set *x);

#endif

