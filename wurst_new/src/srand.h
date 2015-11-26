/*
 * 29 Sep 2004
 * rcsid = $Id: srand.h,v 1.1 2007/05/11 09:54:22 mahmood Exp $
 */

#ifndef SRAND_H
#define SRAND_H

float          f_rand (void);
void           ini_rand (long int seed);
unsigned short s_r_rand (unsigned short min, unsigned short max);
size_t         st_r_rand (size_t min, size_t max);
size_t         st_g_rand (size_t mean, size_t std_dev);
#endif /* SRAND_H */
