/*
 * 27 Aug 2001
 * Routines for declaring, killing, copying and printing matrices.
 *
 * #include <stdlib.h>
 * #include <matrices.h>
 */

#include <float.h>  /* for FLT_MAX */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "e_malloc.h"
#include "matrix.h"
#include "mprintf.h"

#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: matrix.c,v 1.6 2005/04/19 07:12:00 torda Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */


/* ---------------- f_matrix ----------------------------------
 * New two-dimensional matrix in malloc()'d memory
 * It is up to the caller to call kill_f_matrix() to free
 * memory.
 * This version is quite nasty. It returns a matrix full of crap.
 * This is not an error. 
 */
float **
f_matrix (const size_t n_rows, const size_t n_cols)
{
    float **matrix;
    size_t i;
    union {
        int i; float f;
    } crap;
    crap.f = FLT_MAX;
    matrix = E_MALLOC (n_rows * sizeof (matrix[0]));
    matrix [0] = E_MALLOC (n_rows * n_cols * sizeof(matrix[0][0]));
    memset (matrix[0], crap.i, n_rows*n_cols * sizeof (matrix[0][0]));
    for (i = 1; i < n_rows; i++)
        matrix[i] = matrix[i-1] + n_cols;
    return matrix;
}

/* ---------------- copy_f_matrix -----------------------------
 * Copy something allocated by f_matrix().
 */
float **
copy_f_matrix (float **src, const size_t n_rows, const size_t n_cols)
{
    float **dst;
    size_t i;
    dst = E_MALLOC (n_rows * sizeof (dst[0]));
    dst [0] = E_MALLOC (n_rows * n_cols * sizeof(dst[0][0]));
    for (i = 1; i < n_rows; i++)
        dst[i] = dst[i-1] + n_cols;
    memcpy (dst[0], src[0], n_cols * n_rows * sizeof (dst[0][0]));
    return dst;
}

#ifdef want_print_f_matrix
/* ---------------- print_f_matrix ----------------------------
 * Mainly for debugging, so it may be #ifdef'd away.
 */
void
dump_f_matrix (const float **mat, const size_t n_rows, const size_t n_cols)
{
    size_t i, j;
    for (i = 0; i < n_rows; i++) {
        for ( j = 0; j < n_cols; j++)
            mprintf ("%5g ", mat[i][j]);
        mprintf ("\n");
    }
}
#endif /* want_print_f_matrix */

/* ---------------- kill_f_matrix -----------------------------
 */
void
kill_f_matrix ( float **matrix)
{
    if ( ! matrix)
        return;
    free_if_not_null (matrix[0]);
    free (matrix);
}

/* ---------------- uc_matrix ---------------------------------
 * Allocate an unsigned char matrix.
 */
unsigned char **
uc_matrix (const size_t n_rows, const size_t n_cols)
{
    unsigned char **matrix;
    size_t i;
    union {
        int i; float f;
    } crap;
    crap.f = 0.0;  /* Sometimes, for fun, set this to FLT_MAX */
    matrix = E_MALLOC (n_rows * sizeof (matrix[0]));
    matrix [0] = E_MALLOC (n_rows*n_cols*sizeof(matrix[0][0]));
    memset (matrix[0], crap.i, n_rows*n_cols*sizeof (matrix[0][0]));
    for (i = 1; i < n_rows; i++)
        matrix[i] = matrix[i-1] + n_cols;
    return matrix;
}

/* ---------------- kill_uc_matrix ----------------------------
 */
void
kill_uc_matrix ( unsigned char **matrix)
{
    if ( ! matrix)
        return;
    free_if_not_null (matrix[0]);
    free (matrix);
}

#ifdef want_print_uc_matrix
/* ---------------- print_uc_matrix ---------------------------
 * Mainly for debugging, so it may be #ifdef'd away.
 */
void
dump_uc_matrix (unsigned char **mat,
                 const size_t n_rows, const size_t n_cols)
{
    size_t i, j;
    for (i = 0; i < n_rows; i++) {
        for ( j = 0; j < n_cols; j++)
            mprintf ("%5u ", (unsigned)mat[i][j]);
        mprintf ("\n");
    }
}
#endif /* want_print_f_matrix */

/* ---------------- d3_array     ------------------------------
 * Allocate space for a three dimensional array.
 * This will allocate space for an array of anything. The size of
 * the object is passed in via the last argument. It has been
 * tested working on things like an array of structures of funny
 * size.
 * n1, n2 and n3 are the dimensions of the array.
 * The return value should be cast to a pointer of the correct
 * type.
 * We do not know what kind of object we are allocating space
 * for, but the last pointer should be of type char *p, since it
 * seems to be illegal to do pointer arithmetic on void *.
 */
void ***
d3_array( const size_t n1, const size_t n2, const size_t n3, const size_t size)
{
    void ***array;
    void **tmp;
    char *p;
    unsigned i, j;
    const unsigned inc = n3 * (size / sizeof (char));

    array                = E_MALLOC (n1 * sizeof (void **));
    array [0]    = tmp   = E_MALLOC (n1 * n2 * sizeof (void*));
    array [0][0] = p     = E_MALLOC (n1 * n2 * n3 * size);

    for (i = 0; i < n1; i++) {
        array [i] = tmp;
        tmp += n2;
        for (j = 0; j < n2; j++) {
            array[i][j] = p;
            p += inc;
        }
    }

    return array;
}

/* ---------------- kill_3d_array -----------------------------
 * Free the memory associated with a 3 dimensional array.
 */
void
kill_3d_array ( void ***p)
{
    free (p[0][0]);
    free (p[0]);
    free (p);
}
