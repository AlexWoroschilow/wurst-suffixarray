/*
 * 10 January 2002
 * rcsid = $Id: fx.h,v 1.1 2002/02/07 07:33:22 torda Exp $
 */

#ifndef FX_H
#define FX_H

typedef struct FXParam {
    size_t nr_groups, nr_inst, nr_dbins;
    float  *cw;
    float  **Ijk, **Ijk_nbr, **Ijk_psi, *Ijk_dist;
    float  ***paa, ***pna;
    float  **psi, **dpsi, *pdav, *pdsig;
    float  *dbin;
} FXParam;

#endif
