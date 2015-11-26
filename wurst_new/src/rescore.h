/*
 * 5 April 2002
 * rscid = $Id: rescore.h,v 1.2 2002/06/05 07:31:28 torda Exp $
 */
#ifndef RESCORE_H
#define RESCORE_H

struct coord;
float score_rs (struct coord *c, const float *P);

float *param_rs_read (const char *fname);

void  param_rs_destroy (float *p);

#endif
