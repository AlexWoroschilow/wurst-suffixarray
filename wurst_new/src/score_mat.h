/*
 * 12 Sep 2001
 * Define a score matrix. Only people who need to manipulate
 * elements should look here.
 * rcsid = $Id: score_mat.h,v 1.2 2001/10/31 02:33:07 torda Exp $
 */

#ifndef SCORE_MAT_H
#define SCORE_MAT_H

struct score_mat {
    float **mat;
    size_t n_rows, n_cols;
};

#endif /* SCORE_MAT_H */
