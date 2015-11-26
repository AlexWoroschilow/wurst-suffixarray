/*
 * 20 May 2002
 * rcsid = $Id: geo_gap.h,v 1.3 2005/04/29 10:11:53 torda Exp $
 */
#ifndef GEO_GAP_H
#define GEO_GAP_H

struct coord;
int
coord_geo_gap (struct coord *c, float *quad, float *linear, float *logistic,
               unsigned int *num_gap, const float scale, const float max);
#endif /* GEO_GAP_H */
