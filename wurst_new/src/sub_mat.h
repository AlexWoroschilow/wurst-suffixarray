/*
 * 26 March 2001
 * rcsid = $Id: sub_mat.h,v 1.3 2002/03/13 04:52:27 torda Exp $
 */
#ifndef SUB_MAT_H
#define SUB_MAT_H

typedef float mat_t;
struct sub_mat {
    mat_t data [MAX_AA][MAX_AA];
    char *fname;     /* Name of file where we got this from */
    char *comment;   /* Maybe store a comment line or two */
};

#endif /* SUB_MAT_H */
