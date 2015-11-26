/*
 * 27 Jan 2005
 * rcsid = $Id: dipep_i.h,v 1.1 2005/03/02 10:07:59 torda Exp $
 * Warning: this may only be included after <stdlib.h> because
 * we use size_t.
 */
#ifndef DIPEP_H
#define DIPEP_H

struct dpt_list;
struct score_mat;
struct seq;

int              dpt_get_n (const struct dpt_list *dpt_list);
struct dpt_list *dpt_read (const char *fname);
float            dpt_get_val (const struct dpt_list *dpt_list,
                              size_t n, int *error);
void             dpt_list_destroy (struct dpt_list *dpt);
int              dpt_set_val (const struct dpt_list *dpt_list, const size_t n,
                              const float val);
char *           dpt_string (const struct dpt_list *dpt);
int              score_dpt ( struct score_mat *score_mat, struct seq *s1,
                             struct seq *s2, const struct dpt_list *dpt_list);

#endif /* DIPEP_H */
