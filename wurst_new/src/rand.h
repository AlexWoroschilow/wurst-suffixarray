/*
 * 16 nov 2005
 * rcsid = $Id: rand.h,v 1.1 2005/11/17 13:43:36 torda Exp $
 */
#ifndef RAND_H
#define RAND_H

void ini_rand (long int seed);
float g_rand (const float mean, const float std_dev);

#endif /* RAND_H */
