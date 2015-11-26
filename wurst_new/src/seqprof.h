/*
 * 19 nov 2003
 * If you think you have to know about the internal structure of
 * a sequence profile, then you have to include this file.
 * Can only be included after stdlib.h.
 * rcsid = $Id: seqprof.h,v 1.2 2004/02/02 10:25:50 torda Exp $
 */
#ifndef SEQPROF_H
#define SEQPROF_H

enum { blst_afbet_size = 20 };

struct seqprof {
    float **freq_mat;        /* Profile of frequencies */
    struct seq *seq;         /* sequence structure */
    size_t  nres;
};

#endif  /* SEQPROF_H */
