/*
 * 27 Aug 2001
 * rcsid = $Id: matrix.h,v 1.4 2005/04/19 07:12:01 torda Exp $
 */
#ifndef MATRIX_H
#define MATRIX_H

float **f_matrix (const size_t n_rows, const size_t n_cols);
unsigned char **uc_matrix (const size_t n_rows, const size_t n_cols);
void    kill_f_matrix ( float **matrix);
void    kill_uc_matrix (unsigned char **matrix);
float **copy_f_matrix ( float **matrix, const size_t n_rows,
                        const size_t n_cols);

#ifdef  want_print_f_matrix
    void dump_f_matrix (const float **mat, const size_t n_rows,
                         const size_t n_cols);
#endif /* want_print_f_matrix */

#ifdef  want_print_uc_matrix
    void dump_uc_matrix (unsigned char **mat, const size_t n_rows,
                          const size_t n_cols);
#endif /* want_print_uc_matrix */

void ***d3_array( const size_t n1, const size_t n2, const size_t n3,
                  const size_t size);
void kill_3d_array ( void ***p);

#endif /* MATRIX_H */
