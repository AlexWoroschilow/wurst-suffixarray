/*
 * 8 July 2005
 * rcsid = $Id: pair_set_chim.h,v 1.1 2005/07/11 13:13:39 torda Exp $
 */
#ifndef PAIR_SET_CHIM_H
#define PAIR_SET_CHIM_H
struct coord;
struct pair_set;
char *
pair_set_chimera (struct pair_set *pair_set,
                  const struct coord *c1, const struct coord *c2);

#endif /* PAIR_SET_CHIM_H */
